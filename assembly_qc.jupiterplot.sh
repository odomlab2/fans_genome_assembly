#!/bin/bash

#BSUB -q medium
#BSUB -n 4
#BSUB -R rusage[mem=20G]
#BSUB -e assembly_qc.jupiterplot_%I.err
#BSUB -o assembly_qc.jupiterplot_%I.out
#BSUB -J assembly_qc.jupiterplot[1-2]
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate assembly_qc

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
REF=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/data/DMR_v1.0_HiC/Fukomys_damarensis.DMR_v1.0_HiC.fna
ASMS=(
  ${RESULTS_DIR}/scaffolded_assembly/yahs_scaffolds_final.fa
  ${RESULTS_DIR}/curated_assembly/out_JBAT.FINAL.fa
)
ASM=${ASMS[$(( ${LSB_JOBINDEX} - 1 ))]}
REF_PREFIX=$(basename $(dirname ${REF}))
ASM_PREFIX=$(basename $(dirname ${ASM}))
OUT_PREFIX=${ASM_PREFIX}_to_${REF_PREFIX}

mkdir -p ${RESULTS_DIR}/assembly_qc/jupiterplot/${OUT_PREFIX} &&
cd ${RESULTS_DIR}/assembly_qc/jupiterplot/${OUT_PREFIX}

jupiter m=1000000 ng=99 labels=both name=${OUT_PREFIX} ref=${REF} fa=${ASM}