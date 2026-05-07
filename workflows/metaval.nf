/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Flag taxonomy table
include { FLAG_TAXPASTA                                         } from '../modules/local/flag_taxpasta'

// Extract reads of taxIDs
include { TAXID_READS                                           } from '../subworkflows/local/taxid_reads'
include { SEQKIT_FQ2FA as SEQKIT_FQ2FA_READS                    } from '../modules/nf-core/seqkit/fq2fa'
include { PIGZ_UNCOMPRESS                                       } from '../modules/nf-core/pigz/uncompress'
// De novo for extracted taxIDs reads
include { SPADES                                                } from '../modules/nf-core/spades'
include { FLYE                                                  } from '../modules/nf-core/flye'

// BLAST
include { SEQKIT_FQ2FA                                          } from '../modules/nf-core/seqkit/fq2fa'
include { BLAST                                                 } from '../subworkflows/local/blast'
include { BLAST as BLAST_PATHOGEN                               } from '../subworkflows/local/blast'

// Mapping
include { FETCH_BLAST_GENOMES                                   } from '../subworkflows/local/fetch_blast_genomes'
include { MAPPING_SHORTREAD                                     } from '../subworkflows/local/mapping_shortread'
include { MAPPING_LONGREAD                                      } from '../subworkflows/local/mapping_longread'
include { MAPPING_SHORTREAD as MAPPING_SHORTREAD_PATHOGEN       } from '../subworkflows/local/mapping_shortread'
include { MAPPING_LONGREAD as MAPPING_LONGREAD_PATHOGEN         } from '../subworkflows/local/mapping_longread'
include { IGV                                                   } from '../subworkflows/local/igv'
include { IGV as IGV_PATHOGEN                                   } from '../subworkflows/local/igv'

// Metaval static report
include { METAVAL_REPORT                                        } from '../modules/local/metaval_report'

// Calling consensus
include { TAXID_BAM_FASTA as TAXID_BAM_FASTA_SHORTREAD          } from '../subworkflows/local/taxid_bam_fasta'
include { TAXID_BAM_FASTA as TAXID_BAM_FASTA_LONGREAD           } from '../subworkflows/local/taxid_bam_fasta'
include { CONSENSUS                                             } from '../subworkflows/local/consensus'

// Summary subworkflow
include { FASTQC                                                } from '../modules/nf-core/fastqc'
include { MULTIQC                                               } from '../modules/nf-core/multiqc'
include { paramsSummaryMap                                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                } from '../subworkflows/local/utils_nfcore_metaval_pipeline'
include { sample_ntc_branch                                     } from '../subworkflows/local/utils_nfcore_metaval_pipeline'
include { taxpasta_sample_ntc_joined                            } from '../subworkflows/local/utils_nfcore_metaval_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAVAL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_fastqc_files = channel.empty()

    // Filter samplesheet to exclude NTC if params.skip_ntc is true.
    if (params.skip_ntc) {
        ch_samplesheet_filtered = ch_samplesheet.filter { meta, _fastq_1, _fastq_2, _kraken2_report, _kraken2_result, _kraken2_taxpasta, _centrifuge_report, _centrifuge_result, _centrifuge_taxpasta, _diamond, _diamond_taxpasta ->
            return !meta.is_ntc
        }
    } else {
        ch_samplesheet = ch_samplesheet
    }

    // Create input channels for short reads and long reads.
    ch_input = ch_samplesheet_filtered.branch { meta, fastq_1, fastq_2, _kraken2_report, _kraken2_result, _kraken2_taxpasta, _centrifuge_report, _centrifuge_result, _centrifuge_taxpasta, _diamond, _diamond_taxpasta ->

        // Define single_end based on the conditions
        meta.single_end = ( fastq_1 && !fastq_2 )

        // reads channels
        short_reads: meta.instrument_platform != 'OXFORD_NANOPORE'
            return [ meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]

        long_reads: meta.instrument_platform == 'OXFORD_NANOPORE'
            return [ meta, [ fastq_1 ] ]
    }

    //
    // Workflow: Extract reads and verification
    //

    // Verify whether the taxonomic IDs identified by classification are true or false positives.
    if ( params.perform_verify_species ) {

        //
        // LOCAL MODULE: FLAG_TAXPASTA:
        //
        // Flag taxonomy tables by comparing samples with the negative controls (NTC) that have the same meta.library_type and meta.batch.
        if ( params.flag_taxpasta ) {
            // Create taxpasta channels for each classifier
            ch_taxpasta = ch_samplesheet.multiMap { meta, _fastq_1, _fastq_2, _kraken2_report, _kraken2_result, kraken2_taxpasta, _centrifuge_report, _centrifuge_result, centrifuge_taxpasta, _diamond, diamond_taxpasta ->
                kraken2: [ meta + [tool: "kraken2"], kraken2_taxpasta ]
                centrifuge: [ meta + [tool: "centrifuge"], centrifuge_taxpasta ]
                diamond: [ meta + [tool: "diamond"], diamond_taxpasta ]
            }

            // Create sample and NTC taxpasta channels: [ [meta.library_type, meta.batch], meta, taxpasta]
            ch_taxpasta_kraken2 = sample_ntc_branch(ch_taxpasta.kraken2)
            ch_taxpasta_centrifuge =  sample_ntc_branch(ch_taxpasta.centrifuge)
            ch_taxpasta_diamond = sample_ntc_branch(ch_taxpasta.diamond)

            // Join sample and NTC channels for each classifier
            ch_taxpasta_kraken2_ntc = taxpasta_sample_ntc_joined(ch_taxpasta_kraken2.sample, ch_taxpasta_kraken2.ntc)
            ch_taxpasta_centrifuge_ntc = taxpasta_sample_ntc_joined(ch_taxpasta_centrifuge.sample, ch_taxpasta_centrifuge.ntc)
            ch_taxpasta_diamond_ntc = taxpasta_sample_ntc_joined(ch_taxpasta_diamond.sample, ch_taxpasta_diamond.ntc)

            // Create the input channel for FLAG_TAXPASTA
            ch_taxpasta_input = channel.empty()
            ch_taxpasta_input = ch_taxpasta_input.mix(
                ch_taxpasta_kraken2_ntc,
                ch_taxpasta_centrifuge_ntc,
                ch_taxpasta_diamond_ntc
            )
            // The input channels to FLAG_TAXPASTA process: [ meta_sample, taxpasta_sample, meta_ntc, taxpasta_ntc ]
            FLAG_TAXPASTA( ch_taxpasta_input )
        }

        //
        // SUBWORKFLOW: TAXID_READS - extract reads
        //

        // Channels for extracting kraken2/centrifuge/diamond reads
        ch_extract_reads = ch_samplesheet_filtered.multiMap { meta, fastq_1, fastq_2, kraken2_report, kraken2_result, kraken2_taxpasta, centrifuge_report, centrifuge_result, centrifuge_taxpasta, diamond, diamond_taxpasta ->
            kraken2_taxpasta: [ meta + [ tool: "kraken2" ], kraken2_taxpasta ]
            kraken2_report: [ meta + [ tool: "kraken2" ], kraken2_report ]
            kraken2_result: [ meta, kraken2_result ]
            reads:[ meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
            centrifuge_taxpasta: [ meta + [ tool: "centrifuge" ], centrifuge_taxpasta ]
            centrifuge_report: [ meta + [ tool: "centrifuge" ], centrifuge_report ]
            centrifuge_result: [ meta, centrifuge_result ]
            diamond_taxpasta: [ meta + [ tool: "diamond" ], diamond_taxpasta ]
            diamond_tsv: [ meta + [ tool: "diamond" ], diamond ]
        }

        TAXID_READS (
        ch_extract_reads.reads,
        ch_extract_reads.kraken2_taxpasta,
        ch_extract_reads.kraken2_result,
        ch_extract_reads.kraken2_report,
        ch_extract_reads.centrifuge_taxpasta,
        ch_extract_reads.centrifuge_result,
        ch_extract_reads.centrifuge_report,
        ch_extract_reads.diamond_taxpasta,
        ch_extract_reads.diamond_tsv,
        )
        ch_versions            = ch_versions.mix( TAXID_READS.out.versions )

        // Remove empty FASTQ files. This can happen when users want to check if the same species was identified across different classifiers.
        ch_taxid_reads = TAXID_READS.out.reads
            .branch { it ->
                empty: it[0].single_end ? it[1].countFastq() < 1 : it[1][0].countFastq() < 1 || it[1][1].countFastq() < 1
                nonempty: true
            }

        // Transpose to handle both single and paired-end reads
        ch_taxid_reads_transpose = ch_taxid_reads.nonempty
            .map { meta, reads ->
                def read_list = reads instanceof List ? reads: [reads]
                [meta, read_list]
            }
            .transpose ()

        // Convert fastq.gz into fasta files
        SEQKIT_FQ2FA_READS( ch_taxid_reads_transpose )
        PIGZ_UNCOMPRESS ( SEQKIT_FQ2FA_READS.out.fasta )
        //
        // MODULE: DE NOVO - SPADES/FLYE
        //

        // Run de novo assembly if the number of reads exceeds the params.min_read_counts
        ch_taxid_reads_filter = ch_taxid_reads.nonempty
            .branch { meta, reads ->
                blast: meta.single_end ? reads.countFastq() < params.min_read_counts : reads[0].countFastq() < params.min_read_counts || reads[1].countFastq() < params.min_read_counts
                denovo: true
            }
        // Then select first read for BLAST
        ch_blast_reads = ch_taxid_reads_filter.blast
            .map { meta, reads ->
                def read = meta.single_end ? reads : reads[0]
                [ meta, read ]
            }
        // Prepare de novo assembly reads channel for shortreads and longreads
        ch_denovo = ch_taxid_reads_filter.denovo
            .branch { meta, reads ->
                shortreads: meta.instrument_platform != 'OXFORD_NANOPORE'
                    return [ meta, reads, [], [] ]
                longreads: meta.instrument_platform == 'OXFORD_NANOPORE'
                    return [ meta, reads ]
            }
        // Short reads de novo assembly
        ch_contigs_denovo = channel.empty()
        if ( params.perform_shortread_denovo ) {
            SPADES( ch_denovo.shortreads, [], [] )
            ch_contigs_denovo = ch_contigs_denovo.mix( SPADES.out.contigs )
        }
        // Long reads de novo assembly
        if ( params.perform_longread_denovo ) {
            FLYE( ch_denovo.longreads, params.flye_mode )
            ch_contigs_denovo = ch_contigs_denovo.mix( FLYE.out.fasta )
        }

        //
        // SUBWORKFLOW: BLAST
        //

        // Prepare the query fasta file
        if ( (!params.skip_blastn) || (!params.skip_blastx)) {
            SEQKIT_FQ2FA ( ch_blast_reads )
            // Build ch_blast_query fasta file
            ch_blast_query = SEQKIT_FQ2FA.out.fasta
            if ( params.perform_shortread_denovo ) {
                ch_blast_query = ch_blast_query.mix( SPADES.out.contigs )
            }
            if ( params.perform_longread_denovo ) {
                ch_blast_query = ch_blast_query.mix( FLYE.out.fasta )
            }
        }

        BLAST(ch_blast_query, params.blastn_db, params.blastx_db, params.blast_header )

        // Perform FASTQC for reads with BLASTN hits
        ch_fastqc_blastn = ch_taxid_reads.nonempty
            .join( BLAST.out.blastn_filtered, by: 0 )
            .map { meta, reads, _filtered_blast -> [ meta, reads ] }

        ch_fastqc_blastx = ch_taxid_reads.nonempty
            .join( BLAST.out.blastx_filtered, by: 0 )
            .map { meta, reads, _filtered_blast -> [ meta, reads ] }

        ch_fastqc_files = ch_fastqc_files.mix( ch_fastqc_blastn, ch_fastqc_blastx )

        //
        // SUBWORKFLOW: MAPPING
        //

        if (params.perform_mapping) {
            // Fetch genomes of blast hits
            FETCH_BLAST_GENOMES ( params.taxid2genome, BLAST.out.unique_taxid, ch_taxid_reads.nonempty )
            // Mapping - short reads
            ch_mapping_input_shortread = FETCH_BLAST_GENOMES.out.shortreads.join(FETCH_BLAST_GENOMES.out.shortreads_genome, by:0)
            MAPPING_SHORTREAD ( ch_mapping_input_shortread )

            // Mapping - long reads
            ch_mapping_input_longread = FETCH_BLAST_GENOMES.out.longreads.join(FETCH_BLAST_GENOMES.out.longreads_genome, by:0)
            MAPPING_LONGREAD ( ch_mapping_input_longread )

            //
            // SUBWORKFLOW: IGV
            //

            // Filter channels to get bam files which contains mapped reads
            // short reads
            ch_igv_input_shortread = channel.empty()
            ch_igv_input_shortread = ch_igv_input_shortread.mix(MAPPING_SHORTREAD.out.bam)
                .join(MAPPING_SHORTREAD.out.bai, by:0)
                .join(FETCH_BLAST_GENOMES.out.shortreads_genome, by:0)

            // long reads
            ch_igv_input_longread = channel.empty()
            ch_igv_input_longread = ch_igv_input_longread.mix(MAPPING_LONGREAD.out.bam)
                .join(MAPPING_LONGREAD.out.bai, by:0)
                .join(FETCH_BLAST_GENOMES.out.longreads_genome, by:0)

            // Prepare IGV input channels
            ch_igv_input = channel.empty()
            ch_igv_input = ch_igv_input.mix ( ch_igv_input_shortread, ch_igv_input_longread )
            IGV( ch_igv_input )
        }

        //
        // SUBWORKFLOW: STATIC REPORT
        //

        ch_samplesheet_report = channel.fromPath ( params.input, checkIfExists: true )
        // BLASTn results
        ch_blastn_report = channel.empty()
        ch_blastn_report = ch_blastn_report.mix( BLAST.out.blastn_filtered )
        // BLASTx results
        ch_blastx_report = channel.empty()
        ch_blastx_report = ch_blastx_report.mix( BLAST.out.blastx_filtered )
        // Coverage tables
        ch_coverage_tables = channel.empty()
        ch_coverage_tables= ch_coverage_tables.mix( MAPPING_SHORTREAD.out.coverage )
        // Coverager plots
        ch_coverage_plots = channel.empty()
        ch_coverage_plots = ch_coverage_plots.mix( MAPPING_SHORTREAD.out.coverage_plot, MAPPING_LONGREAD.out.coverage_plot )

        // Prepare reads folder for the report, if the reads went through de novo assembly, only include the assembly fasta file in the report.
        ch_reads_fa = channel.empty()
        ch_reads_fa = ch_reads_fa.mix(PIGZ_UNCOMPRESS.out.file)
            .groupTuple(by:0)
            .map { meta, reads -> [ meta, reads ] }

        ch_assembly = channel.empty()
        if ( params.perform_shortread_denovo ) {
            ch_assembly = ch_assembly.mix( SPADES.out.scaffolds.ifEmpty(SPADES.out.contigs) )
        }
        if ( params.perform_longread_denovo ) {
            ch_assembly = ch_assembly.mix( FLYE.out.fasta )
        }

        ch_reads_report = channel.empty()
        ch_reads_report = ch_reads_fa.mix(ch_assembly)
            .groupTuple(by:0)
            .map { meta, files ->
                def assembly = files.find { it instanceof Path }
                def reads = files.find { it instanceof List }
                [ meta, assembly ?: reads ]
            }

        METAVAL_REPORT (
            params.ticket,
            ch_samplesheet_report,
            FLAG_TAXPASTA.out.tsv.collect{it[1]},
            ch_reads_report.collect{it[1]},
            ch_blastn_report.collect{it[1]},
            ch_blastx_report.collect{it[1]},
            ch_coverage_tables.collect{it[1]},
            ch_coverage_plots.collect{it[1]}
        )

    }

    //
    // WORKFLOW: Screen pathogens
    //

    //
    // SUBWORKFLOW: MAPPING
    //
    if ( params.perform_screen_pathogens ) {
        ch_reference = file( params.pathogens_genomes, checkIfExists: true)
        // Map short reads to the pathogens genome
        ch_mapping_pathogen_shortread = ch_input.short_reads
            .map { meta, reads -> [ meta, reads, ch_reference]}
        MAPPING_SHORTREAD_PATHOGEN ( ch_mapping_pathogen_shortread )
        ch_multiqc_files = ch_multiqc_files.mix(MAPPING_SHORTREAD_PATHOGEN.out.mqc)

        // Map long reads to the pathogens genome
        ch_mapping_pathogen_longread = ch_input.long_reads
            .map { meta, reads -> [ meta, reads, ch_reference]}
        MAPPING_LONGREAD_PATHOGEN ( ch_mapping_pathogen_longread )
        ch_multiqc_files = ch_multiqc_files.mix(MAPPING_LONGREAD_PATHOGEN.out.mqc)

        // Subset BAM file for each taxID
        ch_accession2taxid = channel.fromPath ( params.accession2taxid, checkIfExists: true )

        TAXID_BAM_FASTA_SHORTREAD ( MAPPING_SHORTREAD_PATHOGEN.out.bam, MAPPING_SHORTREAD_PATHOGEN.out.bai, ch_accession2taxid, params.min_read_counts )
        TAXID_BAM_FASTA_LONGREAD( MAPPING_LONGREAD_PATHOGEN.out.bam, MAPPING_LONGREAD_PATHOGEN.out.bai, ch_accession2taxid, params.min_read_counts )

        // IGV
        ch_igv_input_pathogen_shortread = TAXID_BAM_FASTA_SHORTREAD.out.taxid_bam
            .join( TAXID_BAM_FASTA_SHORTREAD.out.taxid_bai )
            .map { meta, bam, bai ->
                [ meta, bam, bai, ch_reference ]
            }
        ch_igv_input_pathogen_longread = TAXID_BAM_FASTA_LONGREAD.out.taxid_bam
            .join( TAXID_BAM_FASTA_LONGREAD.out.taxid_bai )
            .map { meta, bam, bai ->
                [ meta, bam, bai, ch_reference ]
            }
        ch_igv_input_pathogen = ch_igv_input_pathogen_shortread.mix( ch_igv_input_pathogen_longread )
        IGV_PATHOGEN ( ch_igv_input_pathogen )

        //
        // SUBWORKFLOW: CONSENSUS - BAM file with the number of mapped reads > params.min_read_counts
        //

        ch_bam_filtered = channel.empty()
        ch_bam_filtered_shortread = TAXID_BAM_FASTA_SHORTREAD.out.taxid_bam.join(TAXID_BAM_FASTA_SHORTREAD.out.taxid_bai, by:0)
        ch_bam_filtered_longread = TAXID_BAM_FASTA_LONGREAD.out.taxid_bam.join(TAXID_BAM_FASTA_LONGREAD.out.taxid_bai, by:0)
        ch_bam_filtered = ch_bam_filtered.mix(ch_bam_filtered_shortread, ch_bam_filtered_longread)

        CONSENSUS ( ch_bam_filtered, [ [], ch_reference ], params.consensus_min_bases )

        // BLAST
        // For pair-end reads, only use read1 for BLAST
        ch_shortread_pathogen_blast_read1 = TAXID_BAM_FASTA_SHORTREAD.out.taxid_fasta
            .filter { _meta, reads ->
                reads[0].countFasta() >= 1 && reads[1].countFasta() >= 1
            }
            .map { meta, reads -> [ meta, reads[0]] }
        ch_longread_pathogen_blast = TAXID_BAM_FASTA_LONGREAD.out.taxid_fasta
            .filter { _meta, reads ->
                reads.countFasta() >= 1
            }
        ch_blast_query_pathogen = ch_shortread_pathogen_blast_read1.mix(
            ch_longread_pathogen_blast,
            CONSENSUS.out.consensus.ifEmpty([])
        )
        BLAST_PATHOGEN ( ch_blast_query_pathogen, params.blastn_db, params.blastx_db, params.blast_header )
    }

    //
    // MODULE: FASTQC
    //
    FASTQC( ch_fastqc_files )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC.out.zip.collect{ zips -> zips[1]}.ifEmpty([]) )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'genomic-medicine-sweden_' + 'metaval_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
