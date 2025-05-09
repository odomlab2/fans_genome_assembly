#!/bin/bash

#BSUB -q verylong
#BSUB -n 40
#BSUB -R rusage[mem=90G]
#BSUB -e gene_annotation.err
#BSUB -o gene_annotation.out
#BSUB -J gene_annotation
#BSUB -N

export DATA_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/data
export RESULTS_DIR=/omics/odcf/analysis/OE0538_projects/DO-0001/fans/results
export SOFTMSK_ASM=${RESULTS_DIR}/repeat_annotation/FukAns_EarlGrey/FukAns_summaryFiles/FukAns.softmasked.fasta
export BRAKER_SIF=${HOME}/programs/braker3.sif

# URLs for proteomes from OrthoDB 12 and NCBI RefSeq
URLS=(
  https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Vertebrata.fa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/064/425/GCF_011064425.1_Rrattus_CSIRO_v1/GCF_011064425.1_Rrattus_CSIRO_v1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/881/025/GCF_016881025.1_HiC_Itri_2/GCF_016881025.1_HiC_Itri_2_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/035/GCF_003704035.1_HU_Pman_2.1.3/GCF_003704035.1_HU_Pman_2.1.3_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/317/375/GCF_000317375.1_MicOch1.0/GCF_000317375.1_MicOch1.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/254/825/GCF_030254825.1_Bangor_MerUng_6.1/GCF_030254825.1_Bangor_MerUng_6.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/664/715/GCF_004664715.2_UCI_PerLeu_2.1/GCF_004664715.2_UCI_PerLeu_2.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/159/225/GCF_023159225.1_ASM2315922v1/GCF_023159225.1_ASM2315922v1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/511/675/GCF_047511675.1_mMarFla1.hap1/GCF_047511675.1_mMarFla1.hap1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/686/445/GCF_902686445.1_mSciCar1.2/GCF_902686445.1_mSciCar1.2_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/762/505/GCF_011762505.1_mArvNil1.pat.X/GCF_011762505.1_mArvNil1.pat.X_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/740/685/GCF_020740685.1_mJacJac1.mat.Y.cur/GCF_020740685.1_mJacJac1.mat.Y.cur_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/947/179/515/GCF_947179515.1_mApoSyl1.1/GCF_947179515.1_mApoSyl1.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/903/995/425/GCF_903995425.1_mOncTor1.1/GCF_903995425.1_mOncTor1.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/095/145/GCF_900095145.1_PAHARI_EIJ_v1.1/GCF_900095145.1_PAHARI_EIJ_v1.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/950/005/125/GCF_950005125.1_mChiNiv1.1/GCF_950005125.1_mChiNiv1.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/903/992/535/GCF_903992535.2_mArvAmp1.2/GCF_903992535.2_mArvAmp1.2_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/903/995/435/GCF_903995435.1_mAcoRus1.1/GCF_903995435.1_mAcoRus1.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/949/786/415/GCF_949786415.1_PerEre_H2_v1/GCF_949786415.1_PerEre_H2_v1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/094/665/GCF_900094665.2_CAROLI_EIJ_v1.1/GCF_900094665.2_CAROLI_EIJ_v1.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/632/895/GCF_008632895.1_UCSF_Mcou_1/GCF_008632895.1_UCSF_Mcou_1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/785/GCF_017639785.1_BCM_Maur_2.0/GCF_017639785.1_BCM_Maur_2.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/190/915/GCF_034190915.1_mCavPor4.1/GCF_034190915.1_mCavPor4.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/247/695/GCF_000247695.1_HetGla_female_1.0/GCF_000247695.1_HetGla_female_1.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/218/885/GCF_021218885.2_Marmota_monax_Labrador192_F-V1.1/GCF_021218885.2_Marmota_monax_Labrador192_F-V1.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/806/735/GCF_902806735.1_Bank_vole1_10x/GCF_902806735.1_Bank_vole1_10x_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/426/925/GCF_003426925.1_ASM342692v1/GCF_003426925.1_ASM342692v1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/622/305/GCF_000622305.1_S.galili_v1.0/GCF_000622305.1_S.galili_v1.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/984/765/GCF_001984765.1_C.can_genome_v1.0/GCF_001984765.1_C.can_genome_v1.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/885/135/GCF_014885135.2_M_Fortis_MF-2015_v1.1/GCF_014885135.2_M_Fortis_MF-2015_v1.1_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/827/085/GCF_007827085.1_ASM782708v3/GCF_007827085.1_ASM782708v3_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/012/274/545/GCF_012274545.1_DMR_v1.0_HiC/GCF_012274545.1_DMR_v1.0_HiC_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/260/255/GCF_000260255.1_OctDeg1.0/GCF_000260255.1_OctDeg1.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/024/711/535/GCF_024711535.1_mDipMer1.0.p/GCF_024711535.1_mDipMer1.0.p_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/907/164/565/GCF_907164565.1_mPsaObe1.curated_primary_1811/GCF_907164565.1_mPsaObe1.curated_primary_1811_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/276/665/GCF_000276665.1_ChiLan1.0/GCF_000276665.1_ChiLan1.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/458/135/GCF_001458135.2_marMar/GCF_001458135.2_marMar_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/785/775/GCF_004785775.1_NIH_TR_1.0/GCF_004785775.1_NIH_TR_1.0_protein.faa.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/885/GCF_000151885.1_Dord_2.0/GCF_000151885.1_Dord_2.0_protein.faa.gz
)

mkdir -p ${DATA_DIR}/protein_sequences && cd ${DATA_DIR}/protein_sequences

# download vertebrate protein sequences generated from OrthoDB v12 and
# annotated protein sequences from several species representing diverse
# Rodent families from NCBI RefSeq
for URL in ${URLS[@]}; do
  wget -q -O - ${URL} | gunzip > $(basename ${URL} .gz)
done

# concatenate protein sequence files
cat *.fa* > proteins.faa

export PROT_SEQ=${DATA_DIR}/protein_sequences/proteins.faa

# create symlinks to STAR generated BAM files with different
# names to avoid braker overwriting them during sorting
mkdir -p ${RESULTS_DIR}/rnaseq_mapping/symlinks
for BAM in ${RESULTS_DIR}/rnaseq_mapping/AS-*-LR-*/*.bam; do
  SAMPLE_NAME=$(basename $(dirname ${BAM}))
  ln -s ${BAM} ${RESULTS_DIR}/rnaseq_mapping/symlinks/${SAMPLE_NAME}.bam
done
# do the same for the BAM indices
for BAM in ${RESULTS_DIR}/rnaseq_mapping/AS-*-LR-*/*.bai; do
  SAMPLE_NAME=$(basename $(dirname ${BAM}))
  ln -s ${BAM} ${RESULTS_DIR}/rnaseq_mapping/symlinks/${SAMPLE_NAME}.bai
done

# get a comma-separated list of BAMS
export BAMS=$(find ${RESULTS_DIR}/rnaseq_mapping/symlinks -name '*.bam' -printf '%p,' | sed 's/,$//')

# copy the Augustus config to a writable location
mkdir -p ${HOME}/programs/augustus
apptainer exec ${BRAKER_SIF} cp -r /opt/Augustus/config ${HOME}/programs/augustus

mkdir -p ${DATA_DIR}/tmp && export APPTAINER_TMPDIR=${DATA_DIR}/tmp/

# run braker3 via apptainer
apptainer exec --env=AUGUSTUS_CONFIG_PATH=${HOME}/programs/augustus/config -B ${DATA_DIR},${RESULTS_DIR} ${BRAKER_SIF} \
  braker.pl \
    --species=Fukomys_anselli \
    --genome=${SOFTMSK_ASM} \
    --bam=${BAMS} \
    --prot_seq=${PROT_SEQ} \
    --busco_lineage=glires_odb10 \
    --workingdir=${RESULTS_DIR}/gene_annotation \
    --gff3 \
    --threads ${LSB_DJOB_NUMPROC}

rm -r ${DATA_DIR}/tmp/