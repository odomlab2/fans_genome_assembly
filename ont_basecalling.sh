#!/bin/bash

#BSUB -q gpu
#BSUB -m gpu-v100-32gb-SXM2
#BSUB -gpu num=4:j_exclusive=yes:gmem=24G
#BSUB -e ont_basecalling.err
#BSUB -o ont_basecalling.out
#BSUB -J ont_basecalling
#BSUB -N

source ~/.bashrc

RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
INPUT_DIRS=(/omics/odcf/project/OE0538/DO-0001/fans/nonOTP/ont/view-by-pid/OE0538_DO-0001_fans_5352/liver01/nanopore_direct_dna/*)

mkdir -p ${RESULTS_DIR}/ont_basecalling

# for each sequencing run
for INPUT_DIR in ${INPUT_DIRS[@]}; do
  
  # create a temporary directory containing symbolic links for all pod5 files (i.e. pass, fail, skip)
  RUN_ID=${INPUT_DIR##*/}
  RUN_DIR=${RESULTS_DIR}/ont_basecalling/pod5_${RUN_ID}
  mkdir -p ${RUN_DIR} && ln -s ${INPUT_DIR}/pod5*/*.pod5 ${RUN_DIR}/
  
  # run basecalling
  dorado basecaller sup ${RUN_DIR}/ --emit-fastq > ${RESULTS_DIR}/ont_basecalling/${RUN_ID}.fastq
  
  # remove temporary directory containing symbolic links
  ( for LINK in ${RUN_DIR}/*; do unlink ${LINK}; done ) && rm -r ${RUN_DIR}

done