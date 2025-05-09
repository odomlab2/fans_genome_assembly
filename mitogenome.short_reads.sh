#!/bin/bash

#BSUB -q long
#BSUB -n 30
#BSUB -R rusage[mem=50G]
#BSUB -e mitogenome.short_reads.err
#BSUB -o mitogenome.short_reads.out
#BSUB -J mitogenome.short_reads
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate mitogenome_short_reads

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
DATA_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/data
SHORT_READS=(/omics/odcf/project/OE0538/DO-0001/fans/sequencing/whole_genome_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/spleen01/paired/run240916_LH00229_0073_A22KJW7LT3/sequence/AS-1334285-LR-72829_R*.fastq.gz)

# mitogenome sequences of Fukomys damarensis (NC_027742.1) and Heterocephalus glaber (NC_015112.1)
SEEDS=${DATA_DIR}/mtDNA_african_mole_rats/african_mole_rats.fasta

# initialize database
get_organelle_config.py -a animal_mt

# assemble mitogenome using mitogenomes of F. damarensis and H. glaber as seeds
get_organelle_from_reads.py \
  -1 ${SHORT_READS[0]} \
  -2 ${SHORT_READS[1]} \
  -s ${SEEDS} \
  --fast \
  -k 21,45,65,85,105 \
  -F animal_mt \
  -o ${RESULTS_DIR}/mitogenome_short_reads \
  -t ${LSB_DJOB_NUMPROC}