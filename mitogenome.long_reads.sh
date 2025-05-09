#!/bin/bash

#BSUB -q long
#BSUB -n 8
#BSUB -R rusage[mem=16G]
#BSUB -e mitogenome.long_reads.err
#BSUB -o mitogenome.long_reads.out
#BSUB -J mitogenome.long_reads
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate mitogenome_long_reads

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
DATA_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/data
LONG_READS=(${RESULTS_DIR}/ont_basecalling/*.fastq)
SHORT_READS=(/omics/odcf/project/OE0538/DO-0001/fans/sequencing/whole_genome_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/spleen01/paired/run240916_LH00229_0073_A22KJW7LT3/sequence/AS-1334285-LR-72829_R*.fastq.gz)
REF=${DATA_DIR}/mtDNA_african_mole_rats/NC_027742.1.fasta

mkdir -p ${RESULTS_DIR}/mitogenome_long_reads && cd ${RESULTS_DIR}/mitogenome_long_reads

# map long reads against reference mitogenome of closely related species
for FASTQ in ${LONG_READS[@]}; do
  minimap2 -x map-ont -t ${LSB_DJOB_NUMPROC} -a -L ${REF} ${FASTQ} | \
    samtools sort -@ ${LSB_DJOB_NUMPROC} - | \
    samtools view -@ ${LSB_DJOB_NUMPROC} -F UNMAP,SECONDARY,SUPPLEMENTARY -b -o $(basename ${FASTQ/.fastq/.ont.sorted.bam}) -
done

# merge long read bams
samtools merge --write-index -@ ${LSB_DJOB_NUMPROC} -o mitoreads.bam##idx##mitoreads.bai *.ont.sorted.bam &&
rm *.ont.sorted.bam

# convert merged bam to fastq
bedtools bamtofastq -i mitoreads.bam -fq mitoreads.fastq

# remove reads with identical IDs (not sure what is creating these duplicated IDs)
seqkit rmdup -d mitoreads.dups.fastq -o mitoreads.uniq.fastq -D mitoreads.dups.txt mitoreads.fastq

# filter reads by length (1-20 kbp)
seqkit seq -m 1000 -M 20000 mitoreads.uniq.fastq > mitoreads.filtered.fastq

# convert fastq to fasta
seqkit fq2fa -o mitoreads.filtered.fasta mitoreads.filtered.fastq

# create local blast database from reads
makeblastdb -in mitoreads.filtered.fasta -dbtype nucl -out mitoreads.filtered.db

# get reference mitogenome size
REF_SIZE=$(grep -v "^>" ${REF} | tr -d '\n' | wc -c)

# set cutoff for percentage of query coverage normalized by subject length
NORM_QCOV_CUTOFF=70

# identify and remove reads that only share a very loose similarity with the reference
blastn -query ${REF} -db mitoreads.filtered.db -outfmt '6 sseqid slen qcovs' | \
  sort -k2,2nr | \
  uniq | \
  awk -v gsize="${REF_SIZE}" '{ printf "%s\t%.2f\n", $0, gsize*$3/$2 }' | \
  awk -v qcov="${NORM_QCOV_CUTOFF}" '$4 > qcov' > mitoreads.blast.out

# get the ids of mitoreads that should be kept
cut -f 1 mitoreads.blast.out > mitoreads.keep

# extract mitoreads filtered by normalized query coverage
seqkit grep -f mitoreads.keep mitoreads.filtered.fastq > mitoreads.filtered2.fastq

# assemble and polish mitogenome from filtered mitoreads
flye \
  --nano-hq mitoreads.filtered2.fastq \
  --threads ${LSB_DJOB_NUMPROC} \
  --iterations 2 \
  --out-dir ${RESULTS_DIR}/mitogenome_long_reads/flye

MT_ASM=${RESULTS_DIR}/mitogenome_long_reads/flye/assembly.fasta

# perform qc on short reads (trim adapters, trim polyG, and
# discard reads shorter than 36 bp and low quality reads)
mkdir -p ${RESULTS_DIR}/mitogenome_long_reads/fastp
SAMPLE=$(basename ${SHORT_READS[0]} _R1.fastq.gz)
fastp \
  --in1 ${SHORT_READS[0]} \
  --in2 ${SHORT_READS[1]} \
  --out1 ${RESULTS_DIR}/mitogenome_long_reads/fastp/$(basename ${SHORT_READS[0]}) \
  --out2 ${RESULTS_DIR}/mitogenome_long_reads/fastp/$(basename ${SHORT_READS[1]}) \
  --detect_adapter_for_pe \
  --length_required 36 \
  --overrepresentation_analysis \
  --json ${RESULTS_DIR}/mitogenome_long_reads/fastp/${SAMPLE}.fastp.json \
  --html ${RESULTS_DIR}/mitogenome_long_reads/fastp/${SAMPLE}.fastp.html \
  --report_title ${SAMPLE} \
  --thread ${LSB_DJOB_NUMPROC}

TRIMMED_READS=(${RESULTS_DIR}/mitogenome_long_reads/fastp/*.fastq.gz)

mkdir -p ${RESULTS_DIR}/mitogenome_long_reads/pilon && cd ${RESULTS_DIR}/mitogenome_long_reads/pilon

# index mitogenome assembly and map trimmed short reads
bwa-mem2 index ${MT_ASM}
bwa-mem2 mem -t ${LSB_DJOB_NUMPROC} -R "@RG\\tID:${SAMPLE}\\tSM:${SAMPLE%-LR-*}" ${MT_ASM} ${TRIMMED_READS[@]} | \
  samtools fixmate -@ ${LSB_DJOB_NUMPROC} -m - - | \
  samtools sort -@ ${LSB_DJOB_NUMPROC} - | \
  samtools markdup -@ ${LSB_DJOB_NUMPROC} - - | \
  samtools view --write-index -@ ${LSB_DJOB_NUMPROC} -F 4 -b -o ${SAMPLE}.bam##idx##${SAMPLE}.bai -

# perform error correction with trimmed short reads
pilon --genome ${MT_ASM} --frags ${SAMPLE}.bam --changes --fix all