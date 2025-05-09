#!/bin/bash

#BSUB -q short
#BSUB -n 1
#BSUB -R rusage[mem=8G]
#BSUB -e contam_screen.filtering.err
#BSUB -o contam_screen.filtering.out
#BSUB -J contam_screen.filtering
#BSUB -N

module load SAMtools

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
ASM=${RESULTS_DIR}/purged_assembly/purged.fa
BUSCO=${RESULTS_DIR}/assembly_qc/compleasm/purged_assembly/glires_odb10/full_table.tsv
BLOBTK_DIR=${RESULTS_DIR}/contamination_screening/blobtoolkit
BLOBTK_TXT=${RESULTS_DIR}/contamination_screening/blobtoolkit.txt

mkdir -p ${RESULTS_DIR}/filtered_assembly && cd ${RESULTS_DIR}/filtered_assembly

# get list of contigs containing busco genes
awk '$2~"Single" || $2~"Duplicated" || $2~"Fragmented" { print $3 }' ${BUSCO} | \
  sort -V | \
  uniq > busco_contigs.txt

# get list of contigs to remove
#
# IMPORTANT:
# if blobtoolkit assigned any contigs as potential contaminants, manually check
# blastn results on 'purged.hits' to discriminate false positives
#
# among the contigs flagged as potential contaminants by blobtoolkit, only the
# ones below do not match either Fukomys damarensis or Heterochephalus glaber
# in the blastn results
cat <<-EOF > contaminants.txt
contig_24890_np1212_1
contig_15910_np1212_1
contig_24891_np1212_1
contig_567_np1212_1
contig_25228_np1212_1
contig_19402_np1212_1
contig_27160_np1212_1
contig_32078_np1212_1
EOF
# conditions: no-hit AND no busco genes
tail -n +2 ${BLOBTK_TXT} | \
  grep 'no-hit' | \
  grep -v -F -f busco_contigs.txt | \
  cut -f 2 > no_hit.no_busco.txt
# conditions: hit AND no busco genes AND length <= 1000 bp
tail -n +2 ${BLOBTK_TXT} | \
  grep -v 'no-hit' | \
  grep -v -F -f busco_contigs.txt | \
  awk '$4<=1000 { print $2 }' > hit.no_busco.1kb_or_shorter.txt
# combine lists
cat contaminants.txt no_hit.no_busco.txt hit.no_busco.1kb_or_shorter.txt | \
  sort -V | \
  uniq > contigs2remove.txt

# filter and index fasta
seqkit grep -n -v -f contigs2remove.txt ${ASM} > filtered.fa
samtools faidx filtered.fa