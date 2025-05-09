#!/bin/bash

#BSUB -q medium
#BSUB -n 1
#BSUB -R rusage[mem=65G]
#BSUB -e arima_hic.index_asm.err
#BSUB -o arima_hic.index_asm.out
#BSUB -J arima_hic.index_asm
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate scaffolding

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/filtered_assembly/filtered.fa

# index genome assembly
bwa-mem2 index ${ASM}