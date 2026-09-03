# genomic-medicine-sweden/metaval: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**genomic-medicine-sweden/metaval** is a bioinformatics pipeline for post-processing of [nf-core/taxprofiler](https://github.com/nf-core/taxprofiler) results and supports two workflows:

- **Verify identified species** extracts reads assigned to viral TaxIDs or a user-defined list of TaxIDs, validates them with BLASTN and/or BLASTX, and can optionally map them to genomes associated with BLAST hits.
- **Pathogen screening** maps reads to a predefined pathogen genome database, separates the mapped reads by pathogen, and validates pathogen-specific reads or consensus sequences with BLAST.

The workflows can be enabled independently or together. Pipeline parameter documentation is generated from `nextflow_schema.json`. Run the following command to see the available parameters:

```bash
nextflow run genomic-medicine-sweden/metaval --help
```

The pipeline, constructed using the `nf-core` [template](https://nf-co.re/tools#creating-a-new-pipeline), utilizing Docker/Singularity containers for easy installation and reproducible results. The implementation follows [Nextflow DSL2](https://docs.seqera.io/nextflow/), employing one container per process for simplified maintenance and dependency management. Processes are sourced from [nf-core/modules](https://github.com/nf-core/modules) for broader accessibility within the Nextflow community.

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Prerequisites

1. Install Nextflow (>=25.10.4) using the instructions [here.](https://nextflow.io/docs/latest/getstarted.html#installation)
2. Install one of the following technologies for full pipeline reproducibility: Docker, Singularity, Podman, Shifter or Charliecloud.

See the [nf-core installation documentation](https://nf-co.re/docs/usage/installation) for setup instructions.

## Input samplesheet

Provide the input CSV with:

```bash
--input /path/to/samplesheet.csv
```

Each row represents one sample. Illumina samples can be single-end or paired-end; Nanopore samples use `fastq_1`.

| Column                | Required                  | Description                                                           |
| --------------------- | ------------------------- | --------------------------------------------------------------------- |
| `sample`              | Yes                       | Unique sample name.                                                   |
| `instrument_platform` | Yes                       | `ILLUMINA` or `OXFORD_NANOPORE`.                                      |
| `library_type`        | Yes                       | `DNA`, `RNA`, or `OTHER`.                                             |
| `is_ntc`              | Yes                       | `true` for a negative control; otherwise `false`.                     |
| `batch`               | Yes                       | Batch identifier used to match samples and negative controls.         |
| `fastq_1`             | Yes                       | Gzipped FASTQ containing read 1, single-end reads, or Nanopore reads. |
| `fastq_2`             | No                        | Gzipped FASTQ containing read 2 for paired-end Illumina data.         |
| `kraken2_report`      | For Kraken2 extraction    | Kraken2 classification report.                                        |
| `kraken2_result`      | For Kraken2 extraction    | Per-read Kraken2 classification output.                               |
| `kraken2_taxpasta`    | For Kraken2 extraction    | Standardized Taxpasta profile.                                        |
| `centrifuge_report`   | For Centrifuge extraction | Centrifuge report in Kraken-style format.                             |
| `centrifuge_result`   | For Centrifuge extraction | Per-read Centrifuge classification output.                            |
| `centrifuge_taxpasta` | For Centrifuge extraction | Standardized Taxpasta profile.                                        |
| `diamond`             | For DIAMOND extraction    | Tab-separated per-read DIAMOND classification output.                 |
| `diamond_taxpasta`    | For DIAMOND extraction    | Standardized Taxpasta profile.                                        |

Example:

```csv title="samplesheet.csv"
sample,instrument_platform,library_type,is_ntc,batch,fastq_1,fastq_2,kraken2_report,kraken2_result,kraken2_taxpasta,centrifuge_report,centrifuge_result,centrifuge_taxpasta,diamond,diamond_taxpasta
sample1,ILLUMINA,DNA,false,batch1,sample1_1.fastq.gz,sample1_2.fastq.gz,sample1.kraken2.kraken2.report.txt,sample1.kraken2.kraken2.classifiedreads.txt,kraken2.tsv,sample1.centrifuge.txt,sample1.centrifuge.results.txt,centrifuge.tsv,sample1.diamond.tsv,diamond.tsv
sample1_ntc,ILLUMINA,DNA,true,batch1,ntc_1.fastq.gz,ntc_2.fastq.gz,ntc.kraken2.kraken2.report.txt,ntc.kraken2.kraken2.classifiedreads.txt,kraken2.tsv,ntc.centrifuge.txt,ntc.centrifuge.results.txt,centrifuge.tsv,ntc.diamond.tsv,diamond.tsv
sample2,OXFORD_NANOPORE,RNA,false,batch2,sample2.fastq.gz,,sample2.kraken2.kraken2.report.txt,sample2.kraken2.kraken2.classifiedreads.txt,kraken2.tsv,sample2.centrifuge.txt,sample2.centrifuge.results.txt,centrifuge.tsv,sample2.diamond.tsv,diamond.tsv
```

### Taxpasta requirements

Taxpasta profiles must contain `taxonomy_id`, `name`, `rank`, and `lineage`, followed by sample abundance columns.

When generating profiles with [nf-core/taxprofiler](https://github.com/nf-core/taxprofiler/blob/main/conf/metaval.config), enable:

```text
--run_profile_standardisation
--taxpasta_add_lineage
--taxpasta_add_rank
--taxpasta_add_name
--taxpasta_taxonomy_dir
```

## Shared reference inputs

### BLASTN database

When BLASTN is enabled, provide a BLAST database directory or `.tar.gz` archive:

```bash
--blastn_db /path/to/blastn_database
```

Databases can be downloaded from the [NCBI BLAST database repository](https://ftp.ncbi.nlm.nih.gov/blast/db/). A virus-focused database can considerably reduce runtime compared with the complete `nt` database.

### BLASTX database

When BLASTX is enabled, supply a DIAMOND database:

```bash
--blastx_db /path/to/database.dmnd
```

See the [DIAMOND documentation](https://github.com/bbuchfink/diamond) for database construction instructions.

### Choosing BLAST modes

BLASTN and BLASTX are enabled by default. Disable either mode with `--skip_blastn` or `--skip_blastx`.

Rules enforced by the pipeline:

- An enabled BLAST mode requires its corresponding database.
- Verify-species mapping requires at least one BLAST mode.
- Pathogen screening requires at least one BLAST mode.
- Both modes may be disabled for verify-species runs when mapping is not enabled, such as when you only want to generate a flagged TAXPASTA table.

## Verify identified species

Enable this workflow with `--perform_verify_species`.

### Selecting TaxIDs automatically

Without `--taxid`, the pipeline extracts detected viral TaxIDs from enabled classifier results. Enable one or more classifiers:

```text
--extract_kraken2_reads
--extract_centrifuge_reads
--extract_diamond_reads
```

Optionally exclude phages, contaminants, or reagent-associated TaxIDs with a one-column file:

```bash
--phages_taxid /path/to/excluded_taxids.txt
```

DIAMOND classifications are filtered using `--evalue_threshold`, which defaults to `0.001`.

### Using user-defined TaxIDs

`--taxid` expects a tab-separated file, not a list of command-line numbers. Each row contains sample name, classifer, taxid and corresponding species name, see example below:

```text
SRR13439799	centrifuge	211044	Influenza A virus (A/Puerto Rico/8/1934(H1N1))
SRR13439799	kraken2	211044	Influenza A virus (A/Puerto Rico/8/1934(H1N1))
SRR13439790	diamond	1920753	Gamaleyavirus
SRR13439813	diamond	1920753	Gamaleyavirus
SRR13439790	kraken2	1920753	Gamaleyavirus
SRR13439813	kraken2	1920753	Gamaleyavirus
SRR13439790	centrifuge	878220	Chryseobacterium sp. StRB126
SRR13439802	centrifuge	878220	Chryseobacterium sp. StRB126
SRR13439813	centrifuge	878220	Chryseobacterium sp. StRB126
SRR13439813	kraken2	878220	Chryseobacterium sp. StRB126
SRR13439813	diamond	878220	Chryseobacterium sp. StRB126
```

Run with:

```bash
--taxid /path/to/taxid_list.tsv
```

When this option is supplied, the user-defined TaxIDs are used instead of automatic viral TaxID selection.

### Negative controls

Taxpasta profiles are compared with negative controls sharing the same `library_type` and `batch`. No separate flagging parameter is required.

`--skip_ntc` defaults to `true`, so negative-control samples are excluded from downstream read extraction and validation. Set it to `false` if the controls should also proceed through downstream analysis:

```bash
--skip_ntc false
```

### De novo assembly

Enable `SPAdes` for Illumina reads with `--perform_shortread_denovo` and enable `Flye` for Nanopore reads with `--perform_longread_denovo`.

Assembly is selected when the extracted-read count reaches `--min_read_counts`, which defaults to `100`. Below the threshold, reads are sent directly to BLAST.

### Mapping BLAST hits

Enable mapping with:

```bash
--perform_mapping
```

Mapping requires a tab-separated `taxid2genome` file:

```bash
--taxid2genome /path/to/taxid2genome.tsv
```

Each row contains:

```text
taxid	organism	path_to_genome
```

For example:

```text
1826872	Candidatus_Nitrosocosmicus_hydrocola_archaea	/path/to/GCA_001870125.1_genomic.fna.gz
2810370	Escherichia_phage_vB_EcoP-ZQ2	/path/to/GCA_019095225.1_genomic.fna.gz
```

The pipeline maps extracted reads to genomes associated with filtered BLAST hits. `Bowtie2` is used for Illumina reads and `minimap2` for Nanopore reads.

### Verify-species example

```bash
nextflow run genomic-medicine-sweden/metaval \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --perform_verify_species \
    --extract_kraken2_reads \
    --extract_centrifuge_reads \
    --extract_diamond_reads \
    --blastn_db /path/to/blastn_db.tar.gz \
    --blastx_db /path/to/diamond.dmnd \
    --perform_shortread_denovo \
    --perform_longread_denovo \
    --perform_mapping \
    --taxid2genome /path/to/taxid2genome.tsv
```

### User-defined TaxID example

```bash
nextflow run genomic-medicine-sweden/metaval \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --perform_verify_species \
    --taxid taxid_list.tsv \
    --extract_kraken2_reads \
    --extract_centrifuge_reads \
    --extract_diamond_reads \
    --blastn_db /path/to/blastn_db.tar.gz \
    --blastx_db /path/to/diamond.dmnd
```

## Pathogen screening

Enable this workflow with:

```bash
--perform_screen_pathogens
```

### Pathogen genome database

Provide a concatenated FASTA containing the pathogen reference genomes:

```bash
--pathogens_genomes /path/to/pathogens.fasta
```

### Accession-to-TaxID map

Provide a tab-separated map associating each FASTA accession with a TaxID and organism name:

```bash
--accession2taxid /path/to/accession2taxid.map
```

Each row contains:

```text
accession	taxid	organism
```

The accession must match the corresponding FASTA sequence identifier.

### Pathogen-specific coverage and depth

Reads are first mapped to the concatenated reference. The resulting BAM is then split by pathogen TaxID. Coverage, depth, and coverage plots are generated separately for each sample and pathogen.

### Consensus calling

Pathogens reaching `--min_read_counts` can be used for consensus calling:

- `--perform_shortread_consensus` enables `samtools consensus` for Illumina reads.
- `--perform_longread_consensus` enables long-read consensus.
- `--longread_consensus_tool medaka` selects Medaka, which is the default.
- `--longread_consensus_tool samtools` selects `samtools consensus`.
- `--consensus_min_bases` sets the minimum retained consensus length and defaults to `50`.

Pathogens below `--min_read_counts` are converted to FASTA and used directly as BLAST input.

### Pathogen-screening example

```bash
nextflow run genomic-medicine-sweden/metaval \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --perform_screen_pathogens \
    --pathogens_genomes /path/to/pathogens.fasta \
    --accession2taxid /path/to/accession2taxid.map \
    --blastn_db /path/to/blastn_db.tar.gz \
    --blastx_db /path/to/diamond.dmnd \
    --perform_shortread_consensus \
    --perform_longread_consensus \
    --longread_consensus_tool medaka
```

## BLAST filtering

Filtered BLAST outputs use the following defaults:

| Parameter             | Default | Meaning                             |
| --------------------- | ------- | ----------------------------------- |
| `--blastn_min_qlen`   | `50`    | Minimum BLASTN query length.        |
| `--blastn_min_pident` | `50`    | Minimum BLASTN percentage identity. |
| `--blastn_min_length` | `50`    | Minimum BLASTN alignment length.    |
| `--blastn_max_evalue` | `0.001` | Maximum BLASTN e-value.             |
| `--blastx_min_qlen`   | `50`    | Minimum BLASTX query length.        |
| `--blastx_min_pident` | `50`    | Minimum BLASTX percentage identity. |
| `--blastx_min_length` | `50`    | Minimum BLASTX alignment length.    |
| `--blastx_max_evalue` | `0.001` | Maximum BLASTX e-value.             |

Metaval uses a fixed tabular BLAST output format for BLASTN and DIAMOND BLASTX so that downstream filtering and the HTML report can parse the same columns. The output columns are defined in `conf/modules.config` and described by the internal header file `assets/blast_outfmt10_header.txt`:

```text
qseqid	staxids	sscinames	pident	qlen	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	sseqid	qseq	sseq
```

Filtering removes hits without a usable `staxids` or `sscinames`, then applies the configured thresholds to `qlen`, `pident`, `length`, and `evalue`.

BLASTX is not run on reads originally extracted from DIAMOND results.

## Parameter files

Pipeline parameters can be supplied using a YAML or JSON file:

```bash
nextflow run genomic-medicine-sweden/metaval \
    -profile docker \
    -params-file params.yaml
```

Example:

```yaml title="params.yaml"
input: samplesheet.csv
outdir: results
perform_verify_species: true
extract_kraken2_reads: true
blastn_db: /path/to/blastn_db.tar.gz
skip_blastx: true
```

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull genomic-medicine-sweden/metaval
```

## Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [genomic-medicine-sweden/metaval releases page](https://github.com/genomic-medicine-sweden/metaval/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend using Docker or Singularity containers for full pipeline reproducibility. When this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please use Conda only as a last resort, i.e., when it is not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from pipeline steps whose inputs are unchanged, continuing from where the previous run stopped. For inputs to be considered identical, both their names and file contents must match. For more information about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Although the pipeline's default requirements should work for most users and input data, you may need to customize the requested compute resources. Each pipeline step has default requirements for CPUs, memory, and execution time. For most steps, if a job exits with one of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18), it is automatically resubmitted with increased resource requests (twice the original request, then three times the original request). If it still fails after the third attempt, pipeline execution stops.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
