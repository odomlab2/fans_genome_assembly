#!/bin/bash

#BSUB -q verylong
#BSUB -n 20
#BSUB -R rusage[mem=100G]
#BSUB -e repeat_annotation.err
#BSUB -o repeat_annotation.out
#BSUB -J repeat_annotation
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate repeat_annotation

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/curated_assembly/out_JBAT.FINAL.fa

mkdir -p ${RESULTS_DIR}/repeat_annotation && cd ${RESULTS_DIR}/repeat_annotation

# annotate repeats based on de novo repeat library and known repeats from Dfam
earlGrey -g ${ASM} -s FukAns -o . -t ${LSB_DJOB_NUMPROC} -r Fukomys -d yes -e yes