#!/bin/bash

#BSUB -q verylong
#BSUB -n 20
#BSUB -R rusage[mem=30G]
#BSUB -e arima_hic.mapping_%I.err
#BSUB -o arima_hic.mapping_%I.out
#BSUB -J arima_hic.mapping[1-4]%2
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate scaffolding

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/filtered_assembly/filtered.fa
TRIMMED_HIC_READS=(${RESULTS_DIR}/scaffolded_assembly/arima_hic/01_fastqs/AS-1333551-LR-*_R*.fastq.gz)
FASTQ=${TRIMMED_HIC_READS[$(( ${LSB_JOBINDEX} - 1 ))]}

mkdir -p ${RESULTS_DIR}/scaffolded_assembly/arima_hic/02_raw_bams &&
cd ${RESULTS_DIR}/scaffolded_assembly/arima_hic/02_raw_bams

# map paired-end reads independently as single-ends
bwa-mem2 mem -t ${LSB_DJOB_NUMPROC} ${ASM} ${FASTQ} | \
  samtools view -@ ${LSB_DJOB_NUMPROC} -b - > $(basename ${FASTQ} .fastq.gz).bam