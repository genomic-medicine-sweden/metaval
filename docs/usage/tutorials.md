# genomic-medicine-sweden/metaval: Tutorials

This page provides guidance on setting up `genomic-medicine-sweden/metaval`. The tutorials walk through the following workflows supported by the pipeline:

- **Verify identified species**, using either automatically detected viral TaxIDs or a user-defined list of TaxIDs.
- **Pathogen screening** against a predefined pathogen genome database.

The pipeline post-processes sequencing reads and classification output from
[nf-core/taxprofiler](https://github.com/nf-core/taxprofiler). It supports
Illumina and Oxford Nanopore reads and classification results from Kraken2,
Centrifuge, and DIAMOND.

## Preparation

### Hardware

The datasets used should be small enough to run on your own laptop or a single server node.

If you wish to use an HPC cluster or cloud environment without running an interactive session through your scheduler, see the [nf-core documentation](https://nf-co.re/docs/usage/configuration#introduction) for guidance on creating an appropriate configuration file.

You will need internet access and at least 1.5 GB of hard-drive space.

### Software

The tutorial assumes that you are using a Unix-based operating system and have already installed Nextflow and a software environment system such as [Conda](https://docs.conda.io/en/latest/miniconda.html), [Docker](https://www.docker.com/), or [Singularity/Apptainer](https://apptainer.org/).

The tutorial uses Docker. However, you can replace references to `docker` with `conda`, `singularity`, or `apptainer`, as appropriate.

### Data

**genomic-medicine-sweden/metaval** is a bioinformatics pipeline for post-processing the results of [nf-core/taxprofiler](https://github.com/nf-core/taxprofiler). In this tutorial, we will use the output files from `nf-core/taxprofiler` for a subset of metagenomic sequencing data, including one Illumina sample and one Nanopore sample. The taxonomy path should be specified when running nf-core/taxprofiler, as this information from the `taxpasta` output will be used by `genomic-medicine-sweden/metaval`. Below is an example of running `nf-core/taxprofiler`. Please check the [usage of nf-core/taxprofiler](https://github.com/nf-core/taxprofiler/tree/main/docs) for detailed instructions on how to run it.

#### Run nf-core/taxprofiler

```bash
nextflow run nf-core/taxprofiler \
  -profile hasta,singularity \
  --input samplesheet.csv \
  --databases databases.csv \
  --outdir taxprofiler_results \
  --perform_shortread_qc \
  --perform_longread_qc \
  --perform_shortread_hostremoval \
  --perform_longread_hostremoval \
  --hostremoval_reference /path/to/host_genome.fna \
  --save_hostremoval_index \
  --save_hostremoval_unmapped \
  --run_kraken2 \
  --kraken2_save_reads \
  --kraken2_save_readclassifications \
  --run_centrifuge \
  --centrifuge_save_reads \
  --run_diamond \
  --run_profile_standardisation \
  --taxpasta_taxonomy_dir /path/to/taxonomy \
  --taxpasta_add_lineage \
  --taxpasta_add_rank \
  --taxpasta_add_name
```

#### Download data

First we will create a directory to run the whole tutorial in.

```bash
mkdir metaval-tutorial
cd metaval-tutorial/

```

Reads could be raw FASTQ files, filtered FASTQ files, or FASTQ files with host genomes removed. In this example, FASTQ files with host genomes removed were used.

```bash
# reads:
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439790_SRR13439790.unmapped_1.fastq.gz
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439790_SRR13439790.unmapped_2.fastq.gz
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439799_SRR13439799.unmapped_other.fastq.gz
# reference
curl -O https://github.com/genomic-medicine-sweden/test-datasets/blob/metaval/reference/reference.fasta.gz
curl -O https://github.com/genomic-medicine-sweden/test-datasets/blob/metaval/reference/accession2taxid.map
# kraken2
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439790_k2_pluspf.kraken2.kraken2.report.txt
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439790_k2_pluspf.kraken2.kraken2.classifiedreads.txt
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439799_k2_pluspf.kraken2.kraken2.report.txt
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439799_k2_pluspf.kraken2.kraken2.classifiedreads.txt
# centrifuge
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439790_p_compressed+h+v.centrifuge.txt
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439790_p_compressed+h+v.centrifuge.results.txt
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439799_p_compressed+h+v.centrifuge.txt
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439799_p_compressed+h+v.centrifuge.results.txt
# diamond
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439790_diamond.diamond.tsv
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/SRR13439799_diamond.diamond.tsv
# taxpasta
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/kraken2_k2_pluspf.tsv
curl -O https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/metaval/testdata/centrifuge_p_compressed+h+v.tsv

```

### Preparing input samplesheet

Provide the nf-core/taxprofiler output files in an input samplesheet in CSV format.
The samplesheet contains 15 columns when using all three supported classifiers: Kraken2, Centrifuge, and DIAMOND.

Create a file named `samplesheet.csv`, copy the following lines into it, and save the file:

```csv title="samplesheet.csv"
sample,instrument_platform,library_type,is_ntc,batch,fastq_1,fastq_2,kraken2_report,kraken2_result,kraken2_taxpasta,centrifuge_report,centrifuge_result,centrifuge_taxpasta,diamond,diamond_taxpasta
SRR13439790,ILLUMINA,DNA,false,batch1,SRR13439790_SRR13439790.unmapped_1.fastq.gz,SRR13439790_SRR13439790.unmapped_2.fastq.gz,SRR13439790_k2_pluspf.kraken2.kraken2.report.txt,SRR13439790_k2_pluspf.kraken2.kraken2.classifiedreads.txt,kraken2_k2_pluspf.tsv,SRR13439790_p_compressed+h+v.centrifuge.txt,SRR13439790_p_compressed+h+v.centrifuge.results.txt,centrifuge_p_compressed+h+v.tsv,SRR13439790_diamond.diamond.tsv,diamond_diamond.tsv
SRR13439799,OXFORD_NANOPORE,OTHER,false,batch3,SRR13439799_SRR13439799.unmapped_other.fastq.gz,,SRR13439799_k2_pluspf.kraken2.kraken2.report.txt,SRR13439799_k2_pluspf.kraken2.kraken2.classifiedreads.txt,kraken2_k2_pluspf.tsv,SRR13439799_p_compressed+h+v.centrifuge.txt,SRR13439799_p_compressed+h+v.centrifuge.results.txt,centrifuge_p_compressed+h+v.tsv,SRR13439799_diamond.diamond.tsv,diamond_diamond.tsv
```

If your nf-core/taxprofiler output files are stored elsewhere, provide their full paths in the corresponding samplesheet columns.

## Running the pipeline

**genomic-medicine-sweden/metaval** supports three tutorial scenarios:

- Pathogen screening
- Verification of automatically identified viral TaxIDs
- Verification of user-defined TaxIDs

The following sections provide an example command for each scenario.

### Tutorial 1: Pathogen screening

Enable this workflow with `--perform_screen_pathogens`. Prepare a pathogen genome reference database and the corresponding accession-to-TaxID map file.

```bash
git clone https://github.com/genomic-medicine-sweden/metaval.git
nextflow run metaval/main.nf \
    -profile singularity \
    --input samplesheet.csv \
    --outdir pathogen_screen_result \
    --perform_screen_pathogens \
    --pathogens_genomes reference.fasta \
    --accession2taxid /path/to/accession2taxid.map \
    --blastn_db /path/to/blastn_database \
    --skip_blastx \
    --perform_longread_consensus \
    --perform_shortread_consensus \
    --longread_consensus_tool 'medaka' \
    --consensus_min_bases 50
```

`--perform_shortread_consensus` uses `samtools consensus` for Illumina reads. For Nanopore reads, `--longread_consensus_tool` accepts `medaka` or `samtools`.

The main results are written to `pathogens/`, including mapping files, pathogen-specific reads, consensus sequences, BLAST results, coverage tables, coverage plots, and IGV reports.

### Tutorial 2: Verify automatically identified viral TaxIDs

Enable `--perform_verify_species` without supplying `--taxid` to identify viral TaxIDs automatically from the selected classifier outputs.

Enable at least one classifier:

```text
--extract_kraken2_reads
--extract_centrifuge_reads
--extract_diamond_reads
```

You can optionally exclude phages, contaminants, or reagent-associated TaxIDs by providing a one-column file:

```text
--phages_taxid /path/to/excluded_taxids.txt
```

Taxpasta profiles are automatically compared with negative controls that share the same `library_type` and `batch`.

By default, `--skip_ntc` is `true`, so negative controls are used for Taxpasta comparison but excluded from downstream extraction and validation. To process the controls downstream as well, use:

```bash
--skip_ntc false
```

The following example enables extraction from all three classifiers, de novo assembly, both BLAST modes, and mapping:

```bash
git clone https://github.com/genomic-medicine-sweden/metaval.git
nextflow run metaval/main.nf \
    -profile singularity \
    --input samplesheet.csv \
    --outdir identified_viruses_results \
    --perform_verify_species \
    --extract_kraken2_reads \
    --extract_centrifuge_reads \
    --extract_diamond_reads \
    --blastn_db /path/to/blastn_database \
    --blastx_db /path/to/diamond_database.dmnd \
    --perform_shortread_denovo \
    --perform_longread_denovo \
    --min_read_counts 20 \
    --perform_mapping \
    --taxid2genome /path/to/taxid2genome.tsv
```

### Tutorial 3: Verify user-defined TaxIDs

Enable this workflow with `--perform_verify_species` and `--taxid`. The user-defined TaxIDs are not limited to viruses and can represent bacteria, fungi, archaea, parasites, or plasmids.

```bash
git clone https://github.com/genomic-medicine-sweden/metaval.git
nextflow run genomic-medicine-sweden/metaval \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --perform_verify_species \
    --taxid /path/to/taxid_list.tsv \
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

## Run with a parameter file

For reproducible runs, pipeline parameters can be stored in a YAML or JSON file and passed with `-params-file`:

```bash
nextflow run genomic-medicine-sweden/metaval \
    -profile docker \
    -params-file params.yaml
```

Example:

```yaml title="params.yaml"
input: samplesheet.csv
outdir: identified_species_results
perform_verify_species: true
extract_kraken2_reads: true
blastn_db: /path/to/blastn_database
skip_blastx: true
```

Use `-params-file` or command-line parameters for pipeline options. The Nextflow `-c` option is intended for execution settings such as resources, executors, containers, and process-specific configuration.

## Resume a run

If a run stops and you correct its inputs or parameters, reuse completed tasks with:

```bash
nextflow run genomic-medicine-sweden/metaval \
    -profile docker \
    -params-file params.yaml \
    -resume
```

Do not change or delete the Nextflow `work` directory before resuming.
