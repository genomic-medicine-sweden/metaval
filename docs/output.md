# genomic-medicine-sweden/metaval: Output

## Introduction

This document describes the output produced by the pipeline. The pipeline contains two main workflows:

- **Verify identified species** extracts reads assigned to viral TaxIDs or a user-defined list of TaxIDs, and validates them with BLASTn and/or BLASTx. It can optionally map the reads to genomes associated with BLAST hits.
- **Pathogen screening** maps reads to a predefined pathogen genome database, separates mapped reads by pathogen, and validates pathogen-specific reads or consensus sequences with BLAST.

The two workflows can be enabled independently or together. Output directories are only created when the corresponding workflow, classifier, or optional analysis step is enabled. All paths below are relative to the top-level results directory.

## Verify identified species

This workflow is enabled with `--perform_verify_species`. It supports two ways of choosing TaxIDs:

- If `--taxid` is not supplied, viral TaxIDs are identified automatically from the supplied classifier results.
- If `--taxid` is supplied, reads are extracted for the user-defined TaxIDs. These TaxIDs are not restricted to viruses.

### Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Decontamination](#decontamination) - Flag taxonomy tables against matched negative controls
- [Extract Viral TaxIDs](#Extract-Viral-TaxIDs) - Extract all viral TaxIDs identified by classifiers.
- [Extract Reads](#Extract-Reads) - Extract reads assigned by Kraken2, Centrifuge, or DIAMOND.
- [De novo assembly](#De-novo-assembly) - Optionally perform de novo assembly.
- [BLAST](#Verify-species-BLAST) - Run BLASTN and/or BLASTX on extracted reads or assemblies.
- [Mapping](#Verify-species-mapping) - Perform mapping against genomes selected from BLAST hits.
- [Coverage and depth](#Verify-species-coverage-and-depth) - Calculate coverage and depth of mapped reads across genomes.
- [IGV reports](#Verify-species-IGV-Reports) - IGV Report for visualizing mapping results.
- [Report](#Static-metaval-Report) - Generate a report summarising the results of the pipeline.
- [FastQC](#fastqc) - QC of reads of TaxIDs assigned by classifiers

### Decontamination

Compare classifier-specific [Taxpasta](https://github.com/taxprofiler/taxpasta) profiles with negative controls that have the same `library_type` and `batch`. This helps distinguish sample-specific classifications from background signal or contamination.

<details markdown="1">
<summary>Output files</summary>

- `taxpasta_flagged/`
  - `<sample_id>_<library_type>_<batch>_<classifier>.tsv`: Flagged taxonomy table for one sample and classifier.

</details>

The comparison column can contain:

- `in_sample`: reads are present only in the sample.
- `in_NTC`: reads are present only in the corresponding negative control.
- `> NTC`: the reads count is higher in the sample than in the corresponding negative control.
- `< NTC`: the reads count is lower in the sample than in the corresponding negative control.
- `equal`: the reads count is equal in the sample and the corresponding negative control.

If no matching negative control is available, the pipeline still creates a sample-specific table without the NTC comparison columns.

The input TAXPASTA table must contain `taxonomy_id`, `name`, `rank`, and `lineage`. When generating the table with [nf-core/taxprofiler](https://github.com/nf-core/taxprofiler), use:

- `--run_profile_standardisation`
- `--taxpasta_add_lineage`
- `--taxpasta_add_rank`
- `--taxpasta_add_name`
- `--taxpasta_taxonomy_dir`

When running [TAXPASTA](https://github.com/taxprofiler/taxpasta) directly, include `--add_lineage`, `--add_rank`, and `--add_name` options.

### Extract Viral TaxIDs

When `--taxid` is not supplied, viral TaxIDs are automatically identified from the selected classifier outputs (Kraken2, Centrifuge, or DIAMOND). TaxIDs listed in `--phages_taxid` are excluded from the results.

<details markdown="1">
<summary>Output files</summary>

- `viral_taxids/`
  - `<sample_id>_kraken2_viral_taxids.tsv`
  - `<sample_id>_centrifuge_viral_taxids.tsv`
  - `<sample_id>_diamond_viral_taxids.tsv`

</details>

Each file contains detected viral TaxIDs and species names. A classifier-specific file is only created when the corresponding extraction option is enabled:

- `--extract_kraken2_reads`
- `--extract_centrifuge_reads`
- `--extract_diamond_reads`

### Extract reads

Retrieve the reads of viral TaxIDs predicted by classifiers or a user-defined list of TaxIDs provided via the `--taxid` option. Enable one or more of:

- `--extract_kraken2_reads`
- `--extract_centrifuge_reads`
- `--extract_diamond_reads`

<details markdown="1">
<summary>Output files</summary>

- `extracted_reads/`
  - `kraken2/`
    - `<sample_id>_taxid_<taxid>_<species>.extracted_kraken2_*.fa`
  - `centrifuge/`
    - `<sample_id>_taxid_<taxid>_<species>.extracted_centrifuge_*.fa`
  - `diamond/`
    - `<sample_id>_taxid_<taxid>_<species>.extracted_diamond_*.fa`

</details>

Only directories for enabled classifiers are created.

### De novo assembly

Extracted reads can be assembled before BLAST:

- `--perform_shortread_denovo` enables `SPAdes` for Illumina reads.
- `--perform_longread_denovo` enables `Flye` for Oxford Nanopore reads.

Assembly is attempted when the number of extracted reads reaches `--min_read_counts`. Below that threshold, extracted reads are used directly as BLAST input.

<details markdown="1">
<summary>Output files</summary>

- `spades/<classifier>/`
  - `<sample_id>_taxid_<taxid>_<species>_<classifier>.contigs.fa`: SPAdes assembly contigs
  - `<sample_id>_taxid_<taxid>_<species>_<classifier>.scaffolds.fa` : SPAdes scaffold assembly
  - `*.log` : SPAdes log or warning log file.

- `flye/<classifier>/`
  - `<sample_id>_taxid_<taxid>_<species>_<classifier>.fasta`: Final assembly in fasta format generated by `Flye`.
  - `*.log`: Log file summarizing steps and intermediate results.

</details>

See the [SPAdes documentation](https://ablab.github.io/spades/) and [Flye documentation](https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md) for detailed output descriptions.

### Verify-species BLAST

BLAST validates the taxonomic assignment of extracted reads or assemblies:

- BLASTN is enabled unless `--skip_blastn` is supplied.
- DIAMOND BLASTX is enabled unless `--skip_blastx` is supplied.

Both modes may be disabled when mapping is not requested. At least one mode must remain enabled when `--perform_mapping` is used because mapping references are selected from BLAST hits.

<details markdown="1">
<summary>Output files</summary>

- `blast/blastn/<classifier>/`
  - `<sample_id>_taxid_<taxid>_<species>_<classifier>_blast_filtered.txt`: Filtered BLASTN hits.
  - `<sample_id>_taxid_<taxid>_<species>_<classifier>_blast_filtered_summary.txt`: Summary of the filtered BLASTN hits.
- `blast/blastx/<classifier>/`
  - `<sample_id>_taxid_<taxid>_<species>_<classifier>_blastx_filtered.txt`: Filtered DIAMOND BLASTX hits.
  - `<sample_id>_taxid_<taxid>_<species>_<classifier>_blastx_filtered_summary.txt`: Summary of the filtered DIAMOND BLASTX hits.

</details>

Provide the relevant databases with `--blastn_db` and `--blastx_db`. Filtering thresholds are controlled by:

- `--blastn_min_qlen`, `--blastn_min_pident`, `--blastn_min_length`, and `--blastn_max_evalue`
- `--blastx_min_qlen`, `--blastx_min_pident`, `--blastx_min_length`, and `--blastx_max_evalue`

Filtered BLASTN and DIAMOND BLASTX tables use the fixed tabular columns configured in `conf/modules.config` and described by `assets/blast_outfmt10_header.txt`: `qseqid`, `staxids`, `sscinames`, `pident`, `qlen`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`, `sseqid`, `qseq`, and `sseq`. Filtering uses `staxids` and `sscinames` to remove incomplete hits, then applies threshold parameters to `qlen`, `pident`, `length`, and `evalue`.

### Verify-species mapping

Mapping is enabled with `--perform_mapping`. Illumina reads are aligned with `Bowtie2` and Nanopore reads with `minimap2`. Reference genomes are selected from the filtered BLAST hits using the map supplied with `--taxid2genome`.

<details markdown="1">
<summary>Output files</summary>

- `mapping/bowtie2/`
  - `build/mappingorganism_<organism>/`: Bowtie2 reference indices.
  - `align/`
    - `<sample_id>_<classifier>_taxid_<taxid>_<species>_mappingorganism_<organism>_<genome_id>_sorted.bam`
    - `<sample_id>_<classifier>_taxid_<taxid>_<species>_mappingorganism_<organism>_<genome_id>_sorted.bam.bai`
- `mapping/minimap2/`
  - `index/mappingorganism_<organism>/`
    - `*.mmi`: minimap2 reference indices.
  - `align/`
    - `<sample_id>_<classifier>_taxid_<taxid>_<species>_mappingorganism_<organism>_<genome_id>_sorted.bam`
    - `<sample_id>_<classifier>_taxid_<taxid>_<species>_mappingorganism_<organism>_<genome_id>_sorted.bam.bai`

</details>

### Verify-species coverage and depth

Coverage and depth outputs are generated when `--perform_mapping` is enabled.

<details markdown="1">
<summary>Output files</summary>

- `samtools/coverage/`
  - `<sample_id>_<classifier>_taxid_<taxid>_<species>_mappingorganism_<organism>_<genome_id>.txt`
- `samtools/depth/`
  - `<sample_id>_<classifier>_taxid_<taxid>_<species>_mappingorganism_<organism>_<genome_id>.tsv`
- `samtools/coverage_plot/`
  - `<sample_id>_<classifier>_taxid_<taxid>_<species>_mappingorganism_<organism>_<genome_id>.png`

</details>

The coverage table summarizes aligned bases and depth across each reference. The depth table contains per-position depth values, and the PNG file visualizes coverage across the genome.

### Verify-species IGV reports

Interactive reports are generated with [igv-reports](https://github.com/igvteam/igv-reports) when mapping is enabled.

<details markdown="1">
<summary>Output files</summary>

- `igv/`
  - `<sample_id>_<classifier>_taxid_<taxid>_<species>_mappingorganism_<organism>_<genome_id>_report.html`: IGV report shows the variants and coverage of reads mapped to the genomes.

</details>

### Static metaval report

The verify-species workflow generates a standalone HTML report. Pathogen-screening results are not currently included.

<details markdown="1">
<summary>Output files</summary>

- `metaval_reports/`
  - `metaval_report.html`: Standalone HTML report for all samples included in the run.

</details>

The report includes:

- Sample metadata from the input samplesheet.
- Classifier-specific Taxpasta tables and NTC comparisons.
- Taxa checked by metaval for each sample and classifier.
- Extracted reads or assemblies used for verification.
- BLASTN and BLASTX result tables.
- Mapping statistics and coverage plots when mapping is enabled.

:::warning
The input file names are used in the `bin/static_report.py` script. If the input file prefixes have changed, please also update them in the `bin/static_reprt.py`
:::

### FastQC

FastQC is run on the extracted verify-species reads that produced filtered BLAST hits when when `--perform_fastqc` is enabled. The default is `false`.

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

## Pathogen screening

The pathogen-screening workflow is enabled with `--perform_screen_pathogens`. Reads are first mapped against the genome database supplied with `--pathogens_genomes`. Mapped reads are then separated by pathogen using `--accession2taxid` and validated with BLAST.

At least one of BLASTN or BLASTX must remain enabled for pathogen screening.

### Pipeline overview - Pathogen screening

The pathogen-screening workflow performs the following steps:

1. [Mapping](#pathogen-screening-mapping) - Map reads to the pathogen genome database.
2. [Coverage and Depth](#pathogen-screening-coverage-and-depth) - Calculate coverage and depth of mapped reads across genomes.
3. [Pathogen specific reads](#pathogen-specific-reads) - Separate mapped reads by pathogen.
4. [Call consensus](#consensus-calling) - Call consensus for pathogens with sufficient mapped reads.
5. [BLAST](#pathogen-screening-blast) - Run BLAST on pathogen reads or consensus sequences.
6. [IGV reports](#pathogen-screening-igv-reports) - IGV Report for visualizing mapping results.

### Pathogen-screening mapping

Illumina reads are aligned with `Bowtie2` and Nanopore reads with `minimap2` against the complete pathogen reference supplied with `--pathogens_genomes`.

<details markdown="1">
<summary>Output files</summary>

- `pathogens/mapping/bowtie2/`
  - `build/`: Bowtie2 reference indices.
  - `align/`
    - `<sample_id>_aligned_pathogens_genome_sorted.bam` : BAM file containing short reads that aligned against the user-provided pathogens genomes
    - `<sample_id>_aligned_pathogens_genome_sorted.bam.bai`: Index of the bam file.
- `pathogens/mapping/minimap2/`
  - `index/`
    - `*.mmi`: minimap2 reference indices.
  - `align/`
    - `<sample_id>_aligned_pathogens_genome_sorted.bam`: BAM file containing long reads that aligned against the user-provided pathogens genomes
    - `<sample_id>_aligned_pathogens_genome_sorted.bam.bai`: Bam file index

</details>

### Pathogen-specific reads

The initial BAM can contain alignments to several pathogen genomes. The pipeline uses `--accession2taxid` to group reference accessions by TaxID and creates pathogen-specific files.

<details markdown="1">
<summary>Output files</summary>

- `pathogens/pathogen_reads/consensus_input/`
  - `<sample_id>_taxid_<taxid>_<species>_sorted.bam`
  - `<sample_id>_taxid_<taxid>_<species>_sorted.bam.bai`
- `pathogens/pathogen_reads/blast_input/`
  - `<sample_id>_taxid_<taxid>_<species>_{1,2}.fasta.gz`

</details>

When the mapped-read count for a pathogen reaches `--min_read_counts`, its BAM file is placed in `consensus_input`. Below that threshold, mapped reads are converted to FASTA and placed in `blast_input`.

### Pathogen-screening coverage and depth

Coverage and depth are separately calculated for each pathogen-specific BAM/reference.

<details markdown="1">
<summary>Output files</summary>

- `pathogens/samtools/coverage/`
  - `<sample_id>_taxid_<taxid>_<species>.txt`: Coverage statistics for each pathogen.
- `pathogens/samtools/depth/`
  - `<sample_id>_taxid_<taxid>_<species>.tsv`: Per-position depth values for each pathogen.
- `pathogens/samtools/coverage_plot/`
  - `<sample_id>_taxid_<taxid>_<species>.png`: Coverage plot for each pathogen.

</details>

### Consensus calling

Consensus calling is optional:

- `--perform_shortread_consensus` uses `samtools consensus` for Illumina reads.
- `--perform_longread_consensus` uses Medaka or `samtools consensus`, selected with `--longread_consensus_tool`.

Consensus sequences shorter than `--consensus_min_bases` are excluded.

<details markdown="1">
<summary>Output files</summary>

- `pathogens/consensus/`
  - `<sample_id>_taxid_<taxid>_<species>_<consensus_tool>.fasta`: Consensus sequence for a pathogen with sufficient mapped reads.

</details>

### Pathogen-screening BLAST

BLAST validates the pathogen-specific reads and consensus sequences:

- Reads below `--min_read_counts` are used directly.
- Consensus sequences can be used for pathogens that reach the read-count threshold.

At least one BLAST mode must remain enabled:

- Supply `--blastn_db` when BLASTN is enabled.
- Supply `--blastx_db` when DIAMOND BLASTX is enabled.

<details markdown="1">
<summary>Output files</summary>

- `pathogens/blast/blastn/`
  - `<sample_id>_taxid_<taxid>_<species>_blast_filtered.txt`: Filtered BLASTN hits.
- `pathogens/blast/blastx/`
  - `<sample_id>_taxid_<taxid>_<species>_blastx_filtered.txt`: Filtered DIAMOND BLASTX hits.

</details>

The same BLAST filtering parameters described for the verify-species workflow are used here.

### Pathogen-screening IGV reports

IGV reports visualize pathogen-specific mapped reads against their references.

<details markdown="1">
<summary>Output files</summary>

- `pathogens/igv/`
  - `<sample_id>_taxid_<taxid>_<species>_report.html`

</details>

## Shared pipeline outputs

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
