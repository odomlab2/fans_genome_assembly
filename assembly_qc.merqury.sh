#!/bin/bash

#BSUB -q long
#BSUB -n 30
#BSUB -R rusage[mem=200G]
#BSUB -e assembly_qc.merqury.err
#BSUB -o assembly_qc.merqury.out
#BSUB -J assembly_qc.merqury
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate assembly_qc

# ensure that the latest versions of 'meryl' and 'merqury' are installed and in your PATH
export PATH=${HOME}/programs/meryl-1.4.1/bin:${PATH}
export MERQURY=${HOME}/programs/merqury

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASMS=(
  ${RESULTS_DIR}/polished_assembly/assembly.fasta
  ${RESULTS_DIR}/corrected_assembly/genome.nextpolish.fasta
  ${RESULTS_DIR}/purged_assembly/purged.fa
  ${RESULTS_DIR}/filtered_assembly/filtered.fa
  ${RESULTS_DIR}/scaffolded_assembly/yahs_scaffolds_final.fa
  ${RESULTS_DIR}/curated_assembly/out_JBAT.FINAL.fa
)
SHORT_READS=(/omics/odcf/project/OE0538/DO-0001/fans/sequencing/whole_genome_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/spleen01/paired/run240916_LH00229_0073_A22KJW7LT3/sequence/AS-1334285-LR-72829_R*.fastq.gz)

mkdir -p ${RESULTS_DIR}/assembly_qc/merqury && cd ${RESULTS_DIR}/assembly_qc/merqury

# build k-mer dbs with meryl
meryl count k=21 threads=${LSB_DJOB_NUMPROC} memory=200 ${SHORT_READS[@]} output read-db.meryl

export OMP_NUM_THREADS=${LSB_DJOB_NUMPROC}

# estimate reference-free quality and k-mer completeness
for ASM in ${ASMS[@]}; do
  ${MERQURY}/merqury.sh read-db.meryl ${ASM} $(basename $(dirname ${ASM}))
done