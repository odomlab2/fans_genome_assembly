[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15489638.svg)](https://doi.org/10.5281/zenodo.15489638)

# Code from: De novo genome assembly of Ansell’s mole-rat (Fukomys anselli)

Code used to process, assemble, and annotate the Ansell's mole-rat genome in an LSF HPC cluster.

## De novo assembly workflow

- Estimation of genome size, heterozygosity, and repetitiveness with [Jellyfish](https://github.com/gmarcais/Jellyfish) and [GenomeScope2.0](https://github.com/tbenavi1/genomescope2.0): `genome_profiling.sh`

- ONT basecalling with [Dorado](https://github.com/nanoporetech/dorado/): `ont_basecalling.sh`

- Quality control of ONT long reads with [NanoPlot](https://github.com/wdecoster/NanoPlot) and [NanoComp](https://github.com/wdecoster/nanocomp): `long_reads_qc.sh`

- De novo assembly and polishing using ONT long reads with [Flye](https://github.com/mikolmogorov/Flye): `polished_assembly.sh`

- SNV/indel error correction using short reads with [NextPolish](https://github.com/Nextomics/NextPolish): `error_correction.sh`

- Removal of haplotigs and overlaps with [purge_dups](https://github.com/dfguan/purge_dups): `purge_dups.sh`

- Contamination screening with [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and [BlobToolKit](https://github.com/genomehubs/blobtoolkit):

  - `contam_screen.mapping.sh`
  - `contam_screen.blastn.sh`
  - `contam_screen.blobtoolkit.sh`
  - `contam_screen.filtering.sh`

- Processing and alignment of Hi-C data to contigs with an adapted version of the [Arima Genomics Hi-C mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline):

  - `arima_hic.fastp.sh`
  - `arima_hic.index_asm.sh`
  - `arima_hic.mapping.sh`
  - `arima_hic.filtering.sh`
  - `arima_hic.pairing.sh`
  - `arima_hic.read_group.sh`
  - `arima_hic.markduplicates.sh`

- Hi-C scaffolding with [YaHS](https://github.com/c-zhou/yahs): `scaffolding.sh`

- Manual curation in [Juicebox Assembly Tools](https://github.com/aidenlab/Juicebox).

- Assessment of assembly contiguity, completeness, and correctness with [QUAST](https://github.com/ablab/quast), [compleasm](https://github.com/huangnengCSU/compleasm), and [Merqury](https://github.com/marbl/merqury):

  - `assembly_qc.quast.sh`
  - `assembly_qc.compleasm.sh`
  - `assembly_qc.merqury.sh`

- Vizualization of assembly synteny with the Damaraland mole-rat assembly ([GCF_012274545.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_012274545.1/)) with [JupiterPlot](https://github.com/JustinChu/JupiterPlot): `assembly_qc.jupiterplot.sh`

## Annotation workflow

### Annotation of repetitive elements

- Fully automated annotation of repeats and transposable elements with [Earl Grey](https://github.com/TobyBaril/EarlGrey): `repeat_annotation.sh`

### Structural gene annotation

- Quality control of RNA-seq data with [fastp](https://github.com/OpenGene/fastp) and [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/):

  - `rnaseq.fastp.sh`
  - `rnaseq.fastq_screen.sh`

- Alignment of RNA-seq data to repeat-softmasked assembly with [STAR](https://github.com/alexdobin/STAR):

  - `rnaseq.index_asm.sh`
  - `rnaseq.mapping.sh`

- Fully automated ab initio and homology-based gene prediction with [BRAKER3](https://github.com/Gaius-Augustus/BRAKER): `gene_annotation.sh`

- Proteome quality assessment with [OMArk](https://github.com/DessimozLab/OMArk): `annotation_qc.omark.sh`

## Mitochondrial genome assembly workflow

- Assemble mitochondrial genome from Illumina PE short reads with [GetOrganelle](https://github.com/Kinggerm/GetOrganelle): `mitogenome.short_reads.sh`

- Assemble mitochondrial genome from ONT long reads: `mitogenome.long_reads.sh`

  1. Extract and filter mitoreads based on the Damaraland mole-rat reference ([NC_027742.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_027742.1)) with [minimap2](https://github.com/lh3/minimap2), [samtools](https://github.com/samtools/samtools), [bedtools](https://bedtools.readthedocs.io/en/latest/index.html), [SeqKit](https://bioinf.shenwei.me/seqkit/), and [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html).
  2. Assemble, polish and error correct Ansell's mole-rat mitogenome with [Flye](https://github.com/mikolmogorov/Flye) and [Pilon](https://github.com/broadinstitute/pilon).

- Adjust starting position and annotate mitogenome with [MitoAnnotator](https://mitofish.aori.u-tokyo.ac.jp/annotation/input/).
