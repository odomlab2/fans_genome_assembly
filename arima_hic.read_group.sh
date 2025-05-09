#!/bin/bash

#BSUB -q long
#BSUB -n 1
#BSUB -R rusage[mem=2G]
#BSUB -e arima_hic.read_group_%I.err
#BSUB -o arima_hic.read_group_%I.out
#BSUB -J arima_hic.read_group[1-2]
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate scaffolding

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
TMP_DIR=${RESULTS_DIR}/scaffolded_assembly/arima_hic/temp
PAIRED_BAMS=(${RESULTS_DIR}/scaffolded_assembly/arima_hic/temp/AS-1333551-LR-*.bam)
PAIRED_BAM=${PAIRED_BAMS[$(( ${LSB_JOBINDEX} - 1 ))]}
SAMPLE=Fans_5352

mkdir -p ${RESULTS_DIR}/scaffolded_assembly/arima_hic/04_paired_bams &&
cd ${RESULTS_DIR}/scaffolded_assembly/arima_hic/04_paired_bams

# add read group
ID=$(basename ${PAIRED_BAM} .bam)
LB=$(echo ${ID} | cut -d '-' -f 2)
picard -Xmx4G -Djava.io.tmpdir=${TMP_DIR} AddOrReplaceReadGroups \
  INPUT=${PAIRED_BAM} \
  OUTPUT=$(basename ${PAIRED_BAM}) \
  ID=${ID} \
  LB=${LB} \
  SM=${SAMPLE} \
  PL=ILLUMINA \
  PU=none