# genomic-medicine-sweden/metaval: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of genomic-medicine-sweden/metaval, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Support for multiple FASTQ files for the same sample. These will be merged at the start of the pipeline []
- Extract taxIDs of viruses
- Extract Kraken2 reads with KrakenTools
- Extract Centrifuge reads
- Extract DIAMOND reads
- de-novo assembly
- BLASTn and BLASTx
- Mapping
- IGV
- Screen pathogens:
  - Map reads to the pathogen genome database using Bowtie2 for short reads and Minimap2 for long reads
  - Call consensus sequences for reads mapped to the pathogen genomes
  - IGV

### `Fixed`

- Removed stale `patch` references for `blast/blastn` and `seqkit/fq2fa` in `modules.json`
- Fixed bowtie2 / minimap2 BAM globs in `tests/default.nf.test` to specific subpaths so each assertion only sees its own outputs
- Removed redundant `[taxid: meta.taxid]` re-merge in `subworkflows/local/taxid_reads/main.nf`

### `Dependencies`

### `Deprecated`
