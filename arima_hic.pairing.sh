#!/bin/bash

#BSUB -q long
#BSUB -n 1
#BSUB -R rusage[mem=4G]
#BSUB -e arima_hic.pairing_%I.err
#BSUB -o arima_hic.pairing_%I.out
#BSUB -J arima_hic.pairing[1-4:2]
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate scaffolding

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
FAIDX=${RESULTS_DIR}/filtered_assembly/filtered.fa.fai
FILT_BAMS=(${RESULTS_DIR}/scaffolded_assembly/arima_hic/03_filtered_bams/AS-1333551-LR-*_R*.bam)
COMBINER=${HOME}/code/mapping_pipeline/two_read_bam_combiner.pl
MAPQ_FILTER=10

mkdir -p ${RESULTS_DIR}/scaffolded_assembly/arima_hic/temp &&
cd ${RESULTS_DIR}/scaffolded_assembly/arima_hic/temp

R1_BAM=${FILT_BAMS[$(( ${LSB_JOBINDEX} - 1 ))]}
R2_BAM=${FILT_BAMS[${LSB_JOBINDEX}]}

# pair, filter, and sort reads
perl ${COMBINER} ${R1_BAM} ${R2_BAM} samtools ${MAPQ_FILTER} | \
  samtools view -b -t ${FAIDX} - | \
  samtools sort -@ ${LSB_DJOB_NUMPROC} -o $(basename ${R1_BAM//_R1}) -