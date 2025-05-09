#!/bin/bash

#BSUB -q long
#BSUB -n 16
#BSUB -R rusage[mem=6G]
#BSUB -o long_reads_qc.out
#BSUB -e long_reads_qc.err
#BSUB -J long_reads_qc
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate long_reads_qc

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
INPUTS=(${RESULTS_DIR}/ont_basecalling/*.fastq)

mkdir -p ${RESULTS_DIR}/long_reads_qc

RUNS=()
# summary statistics and plots for each seqencing run
for INPUT in ${INPUTS[@]}; do
  RUN=$(basename ${INPUT} .fastq)
  RUNS+=(${RUN})
  NanoPlot \
    -t ${LSB_DJOB_NUMPROC} \
    -o ${RESULTS_DIR}/long_reads_qc/${RUN}_nanoplot \
    --loglength \
    --N50 \
    --dpi 300 \
    --fastq ${INPUT}
done

# summary statistics and plots for combined data
NanoPlot \
  -t ${LSB_DJOB_NUMPROC} \
  -o ${RESULTS_DIR}/long_reads_qc/all_nanoplot \
  --loglength \
  --N50 \
  --dpi 300 \
  --fastq ${INPUTS[@]}

# compare different seqencing runs
NanoComp \
  -t ${LSB_DJOB_NUMPROC} \
  -o ${RESULTS_DIR}/long_reads_qc/nanocomp \
  -n ${RUNS[@]} \
  --dpi 300 \
  --fastq ${INPUTS[@]}