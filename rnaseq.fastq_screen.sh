#!/bin/bash

#BSUB -q long
#BSUB -n 16
#BSUB -R rusage[mem=6G]
#BSUB -e rnaseq.fastq_screen.err
#BSUB -o rnaseq.fastq_screen.out
#BSUB -J rnaseq.fastq_screen
#BSUB -N

module load Micromamba

eval "$(micromamba shell hook --shell bash)" && micromamba activate short_reads_qc

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
TRIMMED_READS=(${RESULTS_DIR}/rnaseq_qc/*.fastq.gz)
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/012/274/545/GCF_012274545.1_DMR_v1.0_HiC/GCF_012274545.1_DMR_v1.0_HiC_genomic.fna.gz

cd ${RESULTS_DIR}/rnaseq_qc

# download bowtie2 databases for common genomes
fastq_screen --get_genomes

# download Fukomys damarensis genome
mkdir -p FastQ_Screen_Genomes/F_damarensis && cd FastQ_Screen_Genomes/F_damarensis
wget -q -O - ${URL} | gunzip > Fukomys_damarensis.DMR_v1.0_HiC.fna

# index Fukomys damarensis genome
bowtie2-build Fukomys_damarensis.DMR_v1.0_HiC.fna Fukomys_damarensis.DMR_v1.0_HiC

# manually add genome to 'fastq_screen.conf'

# run contamination screening
fastq_screen \
  --conf ${RESULTS_DIR}/rnaseq_qc/FastQ_Screen_Genomes/fastq_screen.conf \
  --threads ${LSB_DJOB_NUMPROC} \
  ${TRIMMED_READS[@]}