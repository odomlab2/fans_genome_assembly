#!/bin/bash

#BSUB -q long
#BSUB -n 8
#BSUB -R rusage[mem=10G]
#BSUB -e annotation_qc.omark.err
#BSUB -o annotation_qc.omark.out
#BSUB -J annotation_qc.omark
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate annotation_qc

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
PROTEOME=${RESULTS_DIR}/gene_annotation/braker.aa

mkdir -p ${RESULTS_DIR}/annotation_qc/omark && cd ${RESULTS_DIR}/annotation_qc/omark

# download the latest release of the OMA database 
wget https://omabrowser.org/All/LUCA.h5

# install the ete3 NCBI taxonomy database
python -c "from ete3 import NCBITaxa; ncbi = NCBITaxa(); ncbi.update_taxonomy_database()"

# extract and group splicing isoforms of the same gene from the FASTA headers,
# listing all isoforms of each gene on a single line, separated by semicolons.
grep '>' ${PROTEOME} | \
  awk -F'[>.]' '{a[$2] = a[$2] ? a[$2] ";" $2"."$3 : $2"."$3} END {for (i in a) print a[i]}' | \
  sort -V > braker.splice

# generate an OMAmer search result file
omamer search -d LUCA.h5 -q ${PROTEOME} -t ${LSB_DJOB_NUMPROC} -o braker.omamer

# assess proteome completeness, characterize the consistency of all protein
# coding genes with regard to their homologs, and screen for contamination
# based on the OMA orthology database
omark -f braker.omamer -d LUCA.h5 -o . -t 261002 -of ${PROTEOME} -i braker.splice -r family