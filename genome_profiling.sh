#!/bin/bash

#BSUB -q long
#BSUB -n 12
#BSUB -R rusage[mem=180G]
#BSUB -o genome_profiling.out
#BSUB -e genome_profiling.err
#BSUB -J genome_profiling
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate genome_profiling

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
SHORT_READS=(/omics/odcf/project/OE0538/DO-0001/fans/sequencing/whole_genome_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/spleen01/paired/run240916_LH00229_0073_A22KJW7LT3/sequence/AS-1334285-LR-72829_R*.fastq.gz)

mkdir -p ${RESULTS_DIR}/genome_profiling && cd ${RESULTS_DIR}/genome_profiling

# count k-mers
jellyfish count \
  -C \
  -m 21 \
  -s 5G \
  -t ${LSB_DJOB_NUMPROC} \
  -o 21mer.counts \
  -F 2 \
  --quality-start 33 \
  <(zcat ${SHORT_READS[0]}) \
  <(zcat ${SHORT_READS[1]})

# compute histogram of k-mer occurrences
jellyfish histo -t ${LSB_DJOB_NUMPROC} 21mer.counts > 21mer.histo

# estimate genome size, heterozygosity, and repetitiveness
genomescope2 -i 21mer.histo -o . -p 2 -k 21