#!/bin/bash

#BSUB -q long
#BSUB -n 16
#BSUB -R rusage[mem=36G]
#BSUB -e rnaseq.mapping_%I.err
#BSUB -o rnaseq.mapping_%I.out
#BSUB -J rnaseq.mapping[1-18:2]%4
#BSUB -N

module load STAR
module load SAMtools

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
STAR_GENOME_DIR=${RESULTS_DIR}/rnaseq_mapping/STAR_genome
RNASEQ_READS=(${RESULTS_DIR}/rnaseq_qc/*.fastq.gz)
R1=${RNASEQ_READS[$(( ${LSB_JOBINDEX} - 1 ))]}
R2=${RNASEQ_READS[${LSB_JOBINDEX}]}

cd ${RESULTS_DIR}/rnaseq_mapping

# read group
ID=$(basename ${R1} _R1.fastq.gz)

# map RNA-seq reads against repeat-softmasked assembly without annotation
# ('--outSAMstrandField intronMotif' is required for BRAKER)
STAR \
  --runThreadN ${LSB_DJOB_NUMPROC} \
  --genomeDir ${STAR_GENOME_DIR} \
  --readFilesIn ${R1} ${R2} \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:${ID} \
  --outFileNamePrefix ${ID}/ \
  --twopassMode Basic \
  --sjdbOverhang 99 # max(read_length)-1

# index BAM file
samtools index ${RESULTS_DIR}/rnaseq_mapping/${ID}/*.bam