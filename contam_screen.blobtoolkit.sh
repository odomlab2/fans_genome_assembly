#!/bin/bash

#BSUB -q long
#BSUB -n 8
#BSUB -R rusage[mem=10G]
#BSUB -e contam_screen.blobtoolkit.err
#BSUB -o contam_screen.blobtoolkit.out
#BSUB -J contam_screen.blobtoolkit
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate contam_screen

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/purged_assembly/purged.fa
COV=${RESULTS_DIR}/contamination_screening/mapping/Fans.short.bam
HITS=${RESULTS_DIR}/contamination_screening/blastn/purged.hits
TAXDUMP=${RESULTS_DIR}/contamination_screening/taxdump

mkdir -p ${RESULTS_DIR}/contamination_screening && cd ${RESULTS_DIR}/contamination_screening

# create dataset and add coverage and hits
blobtools create \
  --fasta ${ASM} \
  --cov ${COV}=short_read \
  --hits ${HITS} \
  --taxid 261002 \
  --taxdump ${TAXDUMP} \
  --threads ${LSB_DJOB_NUMPROC} \
  blobtoolkit

# set taxonomic rank to order
blobtools replace --key plot.cat=bestsumorder_order blobtoolkit

# create blobplot
blobtools view --plot --format svg blobtoolkit

# create snailplot
blobtools view --plot --view snail --format svg blobtoolkit

# create table with contig info
blobtools filter --table blobtoolkit.txt blobtoolkit

# # if the command above fails, use the one below
# blobtools filter --table STDOUT blobtoolkit > blobtoolkit.json

# # parse 'blobtoolkit.json'using python and pandas
# import json
# import pandas as pd
#
# with open("blobtoolkit.json", "r") as f:
#     data = json.load(f)
# df = pd.DataFrame(data[1:], columns=data[0])
# df.to_csv("blobtoolkit.txt", index=False, sep="\t")