#!/bin/bash

#BSUB -q highmem
#BSUB -n 16
#BSUB -R rusage[mem=300G]
#BSUB -e error_correction.err
#BSUB -o error_correction.out
#BSUB -J error_correction
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate genome_assembly

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/polished_assembly/assembly.fasta
SHORT_READS=(/omics/odcf/project/OE0538/DO-0001/fans/sequencing/whole_genome_sequencing/view-by-pid/OE0538_DO-0001_fans_5352/spleen01/paired/run240916_LH00229_0073_A22KJW7LT3/sequence/AS-1334285-LR-72829_R*.fastq.gz)

mkdir -p ${RESULTS_DIR}/corrected_assembly && cd ${RESULTS_DIR}/corrected_assembly

# prepare file of file names
ls ${SHORT_READS[@]} > sgs.fofn

# create config file
cat <<-EOF > run.cfg
task = best
rewrite = yes
rerun = 3
parallel_jobs = 8
multithread_jobs = ${LSB_DJOB_NUMPROC}
genome = ${ASM}
sgs_fofn = sgs.fofn
workdir = ${RESULTS_DIR}/corrected_assembly
EOF

# run nextpolish
nextPolish run.cfg