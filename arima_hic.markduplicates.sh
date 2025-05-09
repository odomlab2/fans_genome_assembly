#!/bin/bash

#BSUB -q long
#BSUB -n 1
#BSUB -R rusage[mem=40G]
#BSUB -e arima_hic.markduplicates.err
#BSUB -o arima_hic.markduplicates.out
#BSUB -J arima_hic.markduplicates
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate scaffolding

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
PAIR_DIR=${RESULTS_DIR}/scaffolded_assembly/arima_hic/04_paired_bams
TMP_DIR=${RESULTS_DIR}/scaffolded_assembly/arima_hic/temp
STATS=${HOME}/code/mapping_pipeline/get_stats.pl
SAMPLE=Fans_5352

mkdir -p ${RESULTS_DIR}/scaffolded_assembly/arima_hic/05_dedup_bams &&
cd ${RESULTS_DIR}/scaffolded_assembly/arima_hic/05_dedup_bams

# mark duplicates
picard -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=${TMP_DIR} MarkDuplicates \
  $(find ${PAIR_DIR} -name '*.bam' -printf 'INPUT=%p ') \
  OUTPUT=${SAMPLE}.bam \
  METRICS_FILE=${SAMPLE}.metrics.txt \
  TMP_DIR=${TMP_DIR} \
  ASSUME_SORTED=TRUE \
  VALIDATION_STRINGENCY=LENIENT \
  REMOVE_DUPLICATES=TRUE

# index alignment
samtools index ${SAMPLE}.bam

# calculate summary statistics for hi-c links
perl ${STATS} ${SAMPLE}.bam > ${SAMPLE}.bam.stats

rm -r ${TMP_DIR}