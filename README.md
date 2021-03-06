# bovisanalyzer

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**bovisanalyzer** is a bioinformatics best-practice analysis pipeline for nextflow pipeline for analyzing *Mycobacterium bovis* genomes. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Calculate fastq summary statistics ([`fastq-scan`](https://github.com/rpetit3/fastq-scan))
3. Index reference fasta file ([`BWA index`](https://github.com/lh3/bwa))
4. Trim reads for quality and adapter sequence ([`fastp`](https://github.com/OpenGene/fastp))
5. Assign taxonomic labels to sequence reads ([`Kraken 2`](https://ccb.jhu.edu/software/kraken2/))
6. Re-estimate taxonomic abundance of samples analyzed by kraken ([`Bracken`](https://ccb.jhu.edu/software/bracken/))
7. Downsample fastq files ([`Rasusa`](https://github.com/mbhall88/rasusa))
8. Detect drug resistance and lineage ([`TB-Profiler`](https://github.com/jodyphelan/TBProfiler))
9. Detect spoligotype ([`SpoTyping`](https://github.com/xiaeryu/SpoTyping-v2.0))
10. Variant calling
    1. Read mapping ([`BWA mem`](https://github.com/lh3/bwa))
    2. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3. Call and filter variants ([`BCFtools`](http://samtools.github.io/bcftools/bcftools.html))
    4. Convert filtered bcf to pseudogenome fasta ([`vcf2pseudogenome.py`](https://github.com/nf-core/bactmap/blob/dev/bin/vcf2pseudogenome.py))
11. Create alignment from pseudogenome by concatenating fasta files having first checked that the sample sequences are high quality([`calculate_fraction_of_non_GATC_bases.py`](https://github.com/nf-core/bactmap/blob/dev/bin/))
12. Mask alignment using coordinates from [Price *et al* 2018](https://pubmed.ncbi.nlm.nih.gov/30425997/) ([`remove_blocks_from_aln`](https://github.com/sanger-pathogens/remove_blocks_from_aln))
13. Extract variant sites from alignment ([`SNP-sites`](https://github.com/sanger-pathogens/snp-sites))
14. Present QC and visualisation for raw read, alignment, assembly and variant calling results ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download taxonomic databases for Kraken 2 and Bracken (this is a large file and may take a while):

   ```console
   wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz

   tar xvfz minikraken2_v1_8GB_201904.tgz
   ```

4. Create a samplesheet.csv file to use as input for the pipeline. In the directory where your sample fastq files are located, run the following command:

   ```console
   for i in *_1.fastq.gz;do echo ${i%_1.fastq.gz},$(readlink -f $i),$(readlink -f ${i%_1.fastq.gz}_2.fastq.gz) >> samplesheet.csv;done
   ```
   Open `samplesheet.csv` in your favourite terminal text editor (`vim`, `nano`) and add the following line to the top of the file:

   ```console
   sample,fastq_1,fastq_2
   ```

3. Start running your own analysis!

   ```console
   nextflow run avantonder/bovisanalyzer --input samplesheet.csv --reference Mycobacterium_bovis_AF2122_97_GCF_000195835_2.fa -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --kraken2db minikraken2_v1_8GB --brackendb minikraken2_v1_8GB
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

## Documentation

The avantonder/bovisanalyzer pipeline comes with documentation about the pipeline [usage](https://nf-co.re/bovisanalyzer/usage), [parameters](https://nf-co.re/bovisanalyzer/parameters) and [output](https://nf-co.re/bovisanalyzer/output).

## Acknowledgements

avantonder/bovisanalyzer was originally written by Andries van Tonder.

I wouldn't have been able to write this pipeline with out the tools, documentation, pipelines and modules made available by the fantastic [nf-core community](https://nf-co.re/).  In particular bovisanalyzer owes much to the excellent [bactmap](https://nf-co.re/bactmap) pipeline.

## Feedback

If you have any issues, questions or suggestions for improving bovisanalyzer, please submit them to the [Issue Tracker](https://github.com/avantonder/bovisanalyzer/issues).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
