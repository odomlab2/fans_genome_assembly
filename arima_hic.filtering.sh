#!/bin/bash

#BSUB -q long
#BSUB -n 1
#BSUB -R rusage[mem=2G]
#BSUB -e arima_hic.filtering_%I.err
#BSUB -o arima_hic.filtering_%I.out
#BSUB -J arima_hic.filtering[1-4]
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate scaffolding

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
RAW_BAMS=(${RESULTS_DIR}/scaffolded_assembly/arima_hic/02_raw_bams/AS-1333551-LR-*_R*.bam)
RAW_BAM=${RAW_BAMS[$(( ${LSB_JOBINDEX} - 1 ))]}
FILTER=${HOME}/code/mapping_pipeline/filter_five_end.pl

mkdir -p ${RESULTS_DIR}/scaffolded_assembly/arima_hic/03_filtered_bams &&
cd ${RESULTS_DIR}/scaffolded_assembly/arima_hic/03_filtered_bams

# filter 5' end
samtools view -h ${RAW_BAM} | perl ${FILTER} | samtools view -b - > $(basename ${RAW_BAM})