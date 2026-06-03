# genomic-medicine-sweden/metaval

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/genomic-medicine-sweden/metaval)
[![GitHub Actions CI Status](https://github.com/genomic-medicine-sweden/metaval/actions/workflows/nf-test.yml/badge.svg)](https://github.com/genomic-medicine-sweden/metaval/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/genomic-medicine-sweden/metaval/actions/workflows/linting.yml/badge.svg)](https://github.com/genomic-medicine-sweden/metaval/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/genomic-medicine-sweden/metaval)

## Introduction

**genomic-medicine-sweden/metaval** is a bioinformatics pipeline that ...

The pipeline, constructed using the `nf-core` [template](https://nf-co.re/tools#creating-a-new-pipeline), utilizing Docker/Singularity containers for easy installation and reproducible results. The implementation follows [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl1.html), employing one container per process for simplified maintenance and dependency management. Processes are sourced from [nf-core/modules](https://github.com/nf-core/modules) for broader accessibility within the Nextflow community.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,instrument_platform,fastq_1,fastq_2,kraken2_report,kraken2_result,kraken2_taxpasta,centrifuge_report,centrifuge_result,centrifuge_taxpasta,diamond,diamond_taxpasta
sample1,ILLUMINA,sample1.unmapped_1.fastq.gz,sample1.unmapped_2.fastq.gz,sample1.kraken2.kraken2.report.txt,sample1.kraken2.kraken2.classifiedreads.txt,kraken2_kraken2-db.tsv,sample1.centrifuge.txt,sample1.centrifuge.results.txt,centrifuge_centrifuge-db.tsv,sample1.diamond.tsv,diamond_diamond-db.tsv
sample2,ILLUMINA,sample2.unmapped_1.fastq.gz,sample2.unmapped_2.fastq.gz,sample2.kraken2.kraken2.report.txt,sample2.kraken2.kraken2.classifiedreads.txt,kraken2_kraken2-db.tsv,sample2.centrifuge.txt,sample2.centrifuge.results.txt,centrifuge_centrifuge-db.tsv,sample2.diamond.tsv,diamond_diamond-db.tsv
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run genomic-medicine-sweden/metaval \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --perform_verify_species --extract_kraken2_reads
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

genomic-medicine-sweden/metaval was originally written by LilyAnderssonLee.

We thank the following people for their assistance in the development of this pipeline:
[maxulysse](https://github.com/maxulysse)
[mashehu](https://github.com/mashehu)

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use genomic-medicine-sweden/metaval for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
