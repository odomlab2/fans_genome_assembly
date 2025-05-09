#!/bin/bash

#BSUB -q verylong
#BSUB -n 16
#BSUB -R rusage[mem=40G]
#BSUB -e contam_screen.blastn.err
#BSUB -o contam_screen.blastn.out
#BSUB -J contam_screen.blastn
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate contam_screen

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/purged_assembly/purged.fa

mkdir -p ${RESULTS_DIR}/contamination_screening/blastdb && cd ${RESULTS_DIR}/contamination_screening/blastdb

# download the nt and taxbd databases
update_blastdb.pl --source ncbi --decompress --num_threads ${LSB_DJOB_NUMPROC} nt taxdb

mkdir -p ${RESULTS_DIR}/contamination_screening/blastn && cd ${RESULTS_DIR}/contamination_screening/blastn

# perform a blastn search against the nt database
blastn \
  -query ${ASM} \
  -db ${RESULTS_DIR}/contamination_screening/blastdb/nt \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -evalue 1e-25 \
  -outfmt '6 qseqid staxids bitscore std' \
  -num_threads ${LSB_DJOB_NUMPROC} \
  -out $(basename ${ASM} .fa).hits