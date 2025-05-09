#!/bin/bash

#BSUB -q long
#BSUB -n 8
#BSUB -R rusage[mem=12G]
#BSUB -e arima_hic.fastp_%I.err
#BSUB -o arima_hic.fastp_%I.out
#BSUB -J arima_hic.fastp[1-4:2]
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate short_reads_qc

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
HIC_READS=(/omics/odcf/project/OE0538/DO-0001/fans/sequencing/hi_c_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/liver01/paired/run240927_A00372_0870_AHJ3WJDRX5/sequence/AS-1333551-LR-*_R*.fastq.gz)

mkdir -p ${RESULTS_DIR}/scaffolded_assembly/arima_hic/01_fastqs &&
cd ${RESULTS_DIR}/scaffolded_assembly/arima_hic/01_fastqs

R1=${HIC_READS[$(( ${LSB_JOBINDEX} - 1 ))]}
R2=${HIC_READS[${LSB_JOBINDEX}]}
SAMPLE=$(basename ${R1} _R1.fastq.gz)

# trim 5 bases from the 5' end of both reads (skip this if your files are NOT
# prepared with the Arima Hi-C library prep kit), trim adapters, trim polyG,
# and discard reads shorter than 36 bp and low quality reads
fastp \
  --in1 ${R1} \
  --in2 ${R2} \
  --out1 $(basename ${R1}) \
  --out2 $(basename ${R2}) \
  --trim_front1 5 \
  --trim_front2 5 \
  --detect_adapter_for_pe \
  --length_required 36 \
  --overrepresentation_analysis \
  --json ${SAMPLE}.fastp.json \
  --html ${SAMPLE}.fastp.html \
  --report_title ${SAMPLE} \
  --thread ${LSB_DJOB_NUMPROC}