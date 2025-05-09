#!/bin/bash

#BSUB -q long
#BSUB -n 8
#BSUB -R rusage[mem=32G]
#BSUB -e scaffolding.err
#BSUB -o scaffolding.out
#BSUB -J scaffolding
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate scaffolding

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/filtered_assembly/filtered.fa
FAIDX=${ASM}.fai
HIC_BAM=${RESULTS_DIR}/scaffolded_assembly/arima_hic/05_dedup_bams/Fans_5352.bam
JUICER_TOOLS=$HOME/programs/juicer_tools_1.22.01.jar

mkdir -p ${RESULTS_DIR}/scaffolded_assembly && cd ${RESULTS_DIR}/scaffolded_assembly

# run Hi-C scaffolding
yahs -o yahs --no-mem-check ${ASM} ${HIC_BAM}

# generate a Hi-C contact map file for manual curation with Juicebox Assembly Tools (JBAT)
juicer pre -a -o out_JBAT yahs.bin yahs_scaffolds_final.agp ${FAIDX} > out_JBAT.log 2>&1
java -jar -Xmx32G ${JUICER_TOOLS} pre \
  --threads ${LSB_DJOB_NUMPROC} \
  out_JBAT.txt \
  out_JBAT.hic \
  <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{ print $2" "$3 }')

# IMPORTANT:
# Load 'out_JBAT.hic' and 'out_JBAT.assembly' into JBAT for manual editing.
# Once completed editing, export your edits as 'out_JBAT.review.assembly' and
# proceed with the commands below.

# generate AGP and FASTA files for the final genome assembly
juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp ${ASM}

# organize files
mkdir -p ${RESULTS_DIR}/curated_assembly &&
mv out_JBAT.FINAL.{fa,agp} ${RESULTS_DIR}/curated_assembly