#!/bin/bash

#BSUB -q verylong
#BSUB -n 16
#BSUB -R rusage[mem=180G]
#BSUB -e contam_screen.mapping.err
#BSUB -o contam_screen.mapping.out
#BSUB -J contam_screen.mapping
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate contam_screen

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
SHORT_READS=(/omics/odcf/project/OE0538/DO-0001/fans/sequencing/whole_genome_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/spleen01/paired/run240916_LH00229_0073_A22KJW7LT3/sequence/AS-1334285-LR-72829_R*.fastq.gz)
ASM=${RESULTS_DIR}/purged_assembly/purged.fa

mkdir -p ${RESULTS_DIR}/contamination_screening/mapping && cd ${RESULTS_DIR}/contamination_screening/mapping

# map short reads against assembly and mark duplicates
minimap2 -a -x sr -t ${LSB_DJOB_NUMPROC} ${ASM} ${SHORT_READS[0]} ${SHORT_READS[1]} | \
  samtools fixmate -@ ${LSB_DJOB_NUMPROC} -m - - | \
  samtools sort -@ ${LSB_DJOB_NUMPROC} -m 8g - | \
  samtools markdup --write-index -@ ${LSB_DJOB_NUMPROC} - Fans.short.bam