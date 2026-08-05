//
// Subworkflow with functionality specific to the genomic-medicine-sweden/metaval pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
        before_text = """
----------------------------------------------------
   ____ __  __ ____                       _                   _
 / ___|  \\/  / ___|       _ __ ___   ___| |_ __ ___   ____ _| |
| |  _| |\\/| \\___ \\ _____| '_ ` _ \\ / _ \\ __/ _` \\ \\ / / _` | |
| |_| | |  | |___) |_____| | | | | |  __/ || (_| |\\ V / (_| | |
 \\____|_|  |_|____/      |_| |_| |_|\\___| \\__\\__,_| \\_/ \\__,_|_|

  genomic-medicine-sweden/metaval ${workflow.manifest.version}
----------------------------------------------------
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/taxprofiler/blob/main/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    // Fitler NTC or Negative controls from downstream analysis

    ch_samplesheet = channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta,
            fastq_1,
            fastq_2,
            kraken2_report,
            kraken2_result,
            kraken2_taxpasta,
            centrifuge_report,
            centrifuge_result,
            centrifuge_taxpasta,
            diamond,
            diamond_taxpasta ->

            def new_meta = meta + [single_end: fastq_1 && !fastq_2]
            [
            new_meta,
            fastq_1,
            fastq_2,
            kraken2_report,
            kraken2_result,
            kraken2_taxpasta,
            centrifuge_report,
            centrifuge_result,
            centrifuge_taxpasta,
            diamond,
            diamond_taxpasta
            ]
        }

    //
    // Validate parameter inputs
    //

    if (params.perform_verify_species && params.perform_mapping) {
        if (!params.taxid2genome) {
            error ("ERROR: --taxid2genome is required when --perform_mapping is enabled")
        }
        if (params.skip_blastn && params.skip_blastx) {
            error("ERROR: --perform_mapping requires BLASTN or BLASTX. Enable at least one BLAST mode.")
        }
    }

    // At least one BLAST workflow is active.
    def run_blast = params.perform_verify_species || params.perform_screen_pathogens

    if (run_blast && !params.skip_blastn && !params.blastn_db) {
        error ("ERROR: --blastn_db is required when BLASTN is enabled.")
    }

    if (run_blast && !params.skip_blastx && !params.blastx_db) {
        error("ERROR: --blastx_db is required when BLASTX is enabled.")
    }

    if (params.perform_screen_pathogens && params.skip_blastn && params.skip_blastx) {
        error ("ERROR: Pathogen screening requires BLASTN or BLASTX. Enable at least one BLAST mode.")
    }

    if (params.perform_screen_pathogens) {
        if (!params.pathogens_genomes) {
            error ("ERROR: --pathogens_genomes is required with --perform_screen_pathogens.")
        }

        if (!params.accession2taxid) {
            error ("ERROR: --accession2taxid is required with --perform_screen_pathogens.")
        }
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def screen_pathogens = [
        "Mapping reads to a list of viral pathogens genomes with: Bowtie2 (Langmead and Salzberg 2012) for short reads and minimap2 (Li 2018) for long reads,",
        "Consensus calling with either SAMtools (Danecek et al. 2021) or medaka (distributed under the terms of the Oxford Nanopore Technologies PLC. Public License Version 1.0),",
        "Pathogen identification with BLAST (Altschul et al. 1990, Camacho et al. 2009),",
        "Visualisation of results with IGV (Robinson et al. 2011),",
        ].join(' ').trim()

    def verify_species = [
        "Sequencing quality control with: FastQC (Andrews 2010),",
        "Remove false positive findings and background contamination,",
        "Extract reads classified as viruses,",
        "De novo assembly of viral reads with SPAdes (Bankevich et al. 2012) for Illumina reads and Flye (Kolmogorov et al. 2018) for Nanopore reads ,",
        "Reference genome identification with BLAST (Altschul et al. 1990),",
        "Mapping reads to the reference genome with either Bowtie2 (Langmead and Salzberg 2012) for short reads and minimap2 (Li 2018) for long reads,",
        "Depth and coverage statistics with SAMtools (Danecek et al. 2021),",
        "Visualisation of results with IGV (Robinson et al. 2011),",
        ].join(' ').trim()

        def citation_text = [
        "Tools used in the workflow included:",
        params.perform_screen_pathogens ? screen_pathogens : "",
        params.perform_verify_species ? verify_species : "",
        "MultiQC (Ewels et al. 2016)."
    ].join(' ').trim().replaceAll("\\s+", " ").replaceAll("[,|.] +\\.", ".")

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Function that parses and returns the number of mapped reads from flagstat files
//
def getFlagstatMappedReads(flagstat_file) {
    def mapped_reads = 0
    flagstat_file.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped_reads = line.tokenize().first().toInteger()
        }
    }

    def pass = false
   // def logname = flagstat_file.getBaseName() - 'flagstat'
    if (mapped_reads > 0) {
        pass = true
    }
    return [ mapped_reads, pass ]
}

//
// Functions to create input channles for FLAG_TAXPASTA process
//

//Create sample and NTC taxpasta channels
def sample_ntc_branch(taxpasta_channel) {
    return taxpasta_channel.branch { meta, taxpasta ->
        def key = [meta.library_type, meta.batch]
        ntc: meta.is_ntc
        return [key, meta, taxpasta]
        sample: !meta.is_ntc
        return [key, meta, taxpasta]
    }
}

// Join sample and NTC taxpasta channels by meta.library_type and meta.batch
// The input channels to FLAG_TAXPASTA process: [meta_sample, taxpasta_sample, meta_ntc, taxpasta_ntc]
def taxpasta_sample_ntc_joined(ch_taxpasta_sample, ch_taxpasta_ntc) {
    return ch_taxpasta_sample
        .join(ch_taxpasta_ntc, by: 0, remainder: true)
        .filter { entry ->
            // Keep sample-driven rows and discard NTC-only remainder rows from the join
            entry.size() >= 3 && entry[0] != null && entry[1] != null && entry[2] != null
        }
        .map { entry ->
            if (entry.size() == 3) {
                // No matched NTC for this sample; emit sample data with empty NTC placeholders
                def (_key, meta_sample, taxpasta_sample) = entry
                return [meta_sample, taxpasta_sample, [:], []]
            } else {
                // Matched key between sample and NTC rows
                def (_key, meta_sample, taxpasta_sample, meta_ntc, taxpasta_ntc) = entry
                return [meta_sample, taxpasta_sample, meta_ntc ?: [:], taxpasta_ntc ?: []]
            }
        }
}
