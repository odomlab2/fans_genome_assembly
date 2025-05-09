#!/bin/bash

#BSUB -q verylong
#BSUB -n 6
#BSUB -R rusage[mem=20G]
#BSUB -e purge_dups.err
#BSUB -o purge_dups.out
#BSUB -J purge_dups
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate purge_dups

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/corrected_assembly/genome.nextpolish.fasta
LONG_READS=(${RESULTS_DIR}/ont_basecalling/*.fastq)

mkdir -p ${RESULTS_DIR}/purged_assembly && cd ${RESULTS_DIR}/purged_assembly

# run minimap2 to align ont data and generate paf files
for READS in ${LONG_READS[@]}; do
  minimap2 -x map-ont -t ${LSB_DJOB_NUMPROC} -2 ${ASM} ${READS} | \
    pigz -c -p ${LSB_DJOB_NUMPROC} > $(basename ${READS} .fastq).paf.gz
done

# calculate read depth histogram and base-level read depth
pbcstat *.paf.gz
calcuts PB.stat > cutoffs
hist_plot.py -c cutoffs PB.stat PB.cov.png

# split the assembly and do a self-self alignment
SPLIT_ASM=$(basename ${ASM}).split
split_fa ${ASM} > ${SPLIT_ASM}
minimap2 -x asm5 -DP -t ${LSB_DJOB_NUMPROC} -2 ${SPLIT_ASM} ${SPLIT_ASM} | \
  pigz -c -p ${LSB_DJOB_NUMPROC} > ${SPLIT_ASM}.self.paf.gz

# purge haplotigs and overlaps
purge_dups -2 -T cutoffs -c PB.base.cov ${SPLIT_ASM}.self.paf.gz > dups.bed

# get purged primary and haplotig sequences from draft assembly
get_seqs -e dups.bed ${ASM}