#!/bin/bash

#BSUB -q medium
#BSUB -n 16
#BSUB -R rusage[mem=30G]
#BSUB -e rnaseq.index_asm.err
#BSUB -o rnaseq.index_asm.out
#BSUB -J rnaseq.index_asm
#BSUB -N

module load STAR

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
SOFTMSK_ASM=${RESULTS_DIR}/repeat_annotation/FukAns_EarlGrey/FukAns_summaryFiles/FukAns.softmasked.fasta

mkdir -p ${RESULTS_DIR}/rnaseq_mapping && cd ${RESULTS_DIR}/rnaseq_mapping

# index de novo genome assembly without annotation
STAR \
  --runThreadN ${LSB_DJOB_NUMPROC} \
  --runMode genomeGenerate \
  --genomeDir STAR_genome \
  --genomeFastaFiles ${SOFTMSK_ASM}