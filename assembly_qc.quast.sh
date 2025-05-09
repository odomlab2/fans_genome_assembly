#!/bin/bash

#BSUB -q medium
#BSUB -n 12
#BSUB -R rusage[mem=20G]
#BSUB -e assembly_qc.quast.err
#BSUB -o assembly_qc.quast.out
#BSUB -J assembly_qc.quast
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate assembly_qc

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
DATA_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/data
ASMS=(
  ${RESULTS_DIR}/polished_assembly/assembly.fasta
  ${RESULTS_DIR}/corrected_assembly/genome.nextpolish.fasta
  ${RESULTS_DIR}/purged_assembly/purged.fa
  ${RESULTS_DIR}/filtered_assembly/filtered.fa
  ${RESULTS_DIR}/scaffolded_assembly/yahs_scaffolds_final.fa
  ${RESULTS_DIR}/curated_assembly/out_JBAT.FINAL.fa
  ${DATA_DIR}/DMR_v1.0_HiC/Fukomys_damarensis.DMR_v1.0_HiC.fna
  ${DATA_DIR}/HetGla_female_1.0/Heterocephalus_glaber.HetGla_female_1.0.fna
  ${DATA_DIR}/mHetGlaV3/Heterocephalus_glaber.mHetGlaV3.fna
)

mkdir -p ${RESULTS_DIR}/assembly_qc && cd ${RESULTS_DIR}/assembly_qc

# calculate genome assembly statistics
quast \
  -o quast \
  --threads ${LSB_DJOB_NUMPROC} \
  --min-contig 1 \
  --split-scaffolds \
  --large \
  --contig-thresholds 1000,10000,100000,1000000,10000000 \
  -L \
  ${ASMS[@]}