#!/bin/bash

#BSUB -q highmem
#BSUB -n 16
#BSUB -R rusage[mem=850G]
#BSUB -e polished_assembly.err
#BSUB -o polished_assembly.out
#BSUB -J polished_assembly
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate genome_assembly

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
LONG_READS=(${RESULTS_DIR}/ont_basecalling/*.fastq)

flye \
  --nano-hq ${LONG_READS[@]} \
  --genome-size 2.23g \
  --threads ${LSB_DJOB_NUMPROC} \
  --iterations 2 \
  --out-dir ${RESULTS_DIR}/polished_assembly