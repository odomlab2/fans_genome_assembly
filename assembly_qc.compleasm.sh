#!/bin/bash

#BSUB -q long
#BSUB -n 16
#BSUB -R rusage[mem=30G]
#BSUB -e assembly_qc.compleasm.err
#BSUB -o assembly_qc.compleasm.out
#BSUB -J assembly_qc.compleasm
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

mkdir -p ${RESULTS_DIR}/assembly_qc/compleasm && cd ${RESULTS_DIR}/assembly_qc/compleasm

# evaluate genome assemby completeness by testing the presence of a set
# of conserved single-copy orthologs
for ASM in ${ASMS[@]}; do
  compleasm run \
    -a ${ASM} \
    -o $(basename $(dirname ${ASM})) \
    -l glires \
    -t ${LSB_DJOB_NUMPROC}
done