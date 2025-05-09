#!/bin/bash

#BSUB -q medium
#BSUB -n 8
#BSUB -R rusage[mem=4G]
#BSUB -e rnaseq.fastp_%I.err
#BSUB -o rnaseq.fastp_%I.out
#BSUB -J rnaseq.fastp[1-18:2]%4
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate short_reads_qc

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
SHORT_READS=(
  /omics/odcf/project/OE0538/DO-0001/fans/sequencing/rna_sequencing/view-by-pid/OE0538_DO-0001_fans_5342/brain01/paired/run240927_LH00229_0077_B22GCHTLT4/sequence/AS-1333529-LR-72971_R*.fastq.gz
  /omics/odcf/project/OE0538/DO-0001/fans/sequencing/rna_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/lung01/paired/run240927_LH00229_0077_B22GCHTLT4/sequence/AS-1333531-LR-72971_R*.fastq.gz
  /omics/odcf/project/OE0538/DO-0001/fans/sequencing/rna_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/liver01/paired/run240927_LH00229_0077_B22GCHTLT4/sequence/AS-1333533-LR-72971_R*.fastq.gz
  /omics/odcf/project/OE0538/DO-0001/fans/sequencing/rna_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/spleen01/paired/run240927_LH00229_0077_B22GCHTLT4/sequence/AS-1333535-LR-72971_R*.fastq.gz
  /omics/odcf/project/OE0538/DO-0001/fans/sequencing/rna_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/testis01/paired/run240927_LH00229_0077_B22GCHTLT4/sequence/AS-1333537-LR-72971_R*.fastq.gz
  /omics/odcf/project/OE0538/DO-0001/fans/sequencing/rna_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/muscle01/paired/run240927_LH00229_0077_B22GCHTLT4/sequence/AS-1333539-LR-72971_R*.fastq.gz
  /omics/odcf/project/OE0538/DO-0001/fans/sequencing/rna_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/kidney01/paired/run240927_LH00229_0077_B22GCHTLT4/sequence/AS-1333541-LR-72971_R*.fastq.gz
  /omics/odcf/project/OE0538/DO-0001/fans/sequencing/rna_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/heart01/paired/run240927_LH00229_0077_B22GCHTLT4/sequence/AS-1333543-LR-72971_R*.fastq.gz
  /omics/odcf/analysis/OE0538_projects/DO-0001/fans/data/back_skin_RNA/AS-872156-LR-66131_R*.fastq.gz
)

mkdir -p ${RESULTS_DIR}/rnaseq_qc && cd ${RESULTS_DIR}/rnaseq_qc

R1=${SHORT_READS[$(( ${LSB_JOBINDEX} - 1 ))]}
R2=${SHORT_READS[${LSB_JOBINDEX}]}
SAMPLE=$(basename ${R1} _R1.fastq.gz)

# trim adapters, trim polyG, and discard reads shorter than 36 bp and low quality reads
fastp \
  --in1 ${R1} \
  --in2 ${R2} \
  --out1 $(basename ${R1}) \
  --out2 $(basename ${R2}) \
  --detect_adapter_for_pe \
  --length_required 36 \
  --overrepresentation_analysis \
  --json ${SAMPLE}.fastp.json \
  --html ${SAMPLE}.fastp.html \
  --report_title ${SAMPLE} \
  --thread ${LSB_DJOB_NUMPROC}