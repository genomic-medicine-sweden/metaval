//
// Prepare an individual BAM/FASTA file for each pathogen with mapped reads
//

include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PASS               } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FAIL               } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_PASS               } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FAIL               } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX as  SAMTOOLS_INDEX_PASS            } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as  SAMTOOLS_INDEX_FAIL            } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_IDXSTATS                                 } from '../../../modules/nf-core/samtools/idxstats'
include { SAMTOOLS_FASTA                                    } from '../../../modules/nf-core/samtools/fasta'
include { SAMTOOLS_FLAGSTAT                                 } from '../../../modules/nf-core/samtools/flagstat'
include { getFlagstatMappedReads                            } from '../utils_nfcore_metaval_pipeline'
include { SAMTOOLS_COVERAGE                                 } from '../../../modules/nf-core/samtools/coverage'
include { SAMTOOLS_DEPTH                                    } from '../../../modules/nf-core/samtools/depth'
include { COVERAGE_PLOT                                     } from '../../../modules/local/coverage_plot'

workflow TAXID_BAM_FASTA {
    take:
    ch_bam               // Channel: [ val(meta), path(bam) ]
    ch_bai               // Channel: [ val(meta), path(bai) ]
    accession2taxid      // Channel: path(accession2taxid)
    min_read_counts      // Value: minimum number of reads to keep a BAM file

    main:
    ch_taxid_bam       = channel.empty()
    ch_taxid_bai       = channel.empty()
    ch_consensus_input = channel.empty()
    ch_blast_input     = channel.empty()

    // Combine BAM and BAI files
    input_bam = ch_bam.join( ch_bai, by: 0 )
    // Get idxstats for input BAM
    SAMTOOLS_IDXSTATS( input_bam )
    // Extract accessions with mapped reads
    ch_accession_with_meta = SAMTOOLS_IDXSTATS.out.idxstats
        .flatMap { meta, idxstats ->
            idxstats.splitCsv( header: false, sep: "\t" )
                // The SAMTOOLS_IDXSTATS.out.idxstats file contains four columns: <reference_name> <reference_length> <mapped_reads> <unmapped_reads>
                // The last row is "* 0 0 0" and should be filtered out, along with rows that have zero mapped reads.
                .findAll { row -> row[0] != "*" && row[2].toInteger() > 0 }
                .collect{ row -> [meta, row[0], row[2].toInteger()] }
        }

    // Load accession2taxid map
    ch_accession2taxidmap = accession2taxid.splitCsv( header: false, sep: "\t" )

    // Join accessions with taxids: [meta, accession, num_reads] + [accession, taxid, organism]
    ch_accession_taxid_with_meta = ch_accession_with_meta
        .combine(ch_accession2taxidmap)
        .filter { _meta, accession, _num_reads, ref_accession, _taxid, _organism ->
            accession == ref_accession
        }
        .map { meta, accession, num_reads, _ref_accession, taxid, organism ->
            [meta, accession, taxid, organism, num_reads]
        }
        .groupTuple( by: [0, 2, 3] ) // Group by [meta, taxid, organism]
        .map { meta, accession_list, taxid, organism, num_reads_list ->
            [meta, accession_list, taxid, organism, num_reads_list.sum()]
        }
        .branch { meta, accession_list, taxid, organism, num_reads_list ->
            pass: num_reads_list >= min_read_counts // The number of mapped reads to a taxID greater than params.min_read_counts
                return [meta, accession_list, taxid, organism] // [meta, accession_list, taxid, organism]
            fail: num_reads_list < min_read_counts  // The number of mapped reads to a taxID smaller than params.min_read_counts
                return [meta, accession_list, taxid, organism] // [meta, accession_list, taxid, organism]
        }

    // Prepare individual BAM files for each taxID with the number of mapped reads greater than params.min_read_counts
    ch_consensus_input = ch_accession_taxid_with_meta.pass
        .join( input_bam, by: 0 ) // Join by meta (index 0)
        .map { meta, accession_list, taxid, organism, bam, bam_index ->
            // Create new meta with taxid and organism information
            def new_meta = meta + [taxid: taxid, organism: organism, accessions: accession_list.flatten()]
            return [ new_meta, bam, bam_index ]
        }

    // BAM files will be used to call consensus sequences
    SAMTOOLS_VIEW_PASS(ch_consensus_input, [[],[],[]], [[],[]], [[],[]], [] )
    
    // Drop the transient 'accessions' key so downstream meta matches output
    ch_pass_subset = SAMTOOLS_VIEW_PASS.out.bam
        .map { meta, bam -> [meta.subMap(meta.keySet() - 'accessions'), bam ] }
    
    SAMTOOLS_SORT_PASS( ch_pass_subset, [[],[],[]], 'bai' )
    SAMTOOLS_INDEX_PASS( SAMTOOLS_SORT_PASS.out.bam )

    // samtools flagstat check if there are any reads mapped to the genome
    SAMTOOLS_FLAGSTAT (SAMTOOLS_SORT_PASS.out.bam.join(SAMTOOLS_INDEX_PASS.out.index))

    ch_mapped_reads = SAMTOOLS_FLAGSTAT.out.flagstat
        .map { meta, flagstat -> [meta] + getFlagstatMappedReads(flagstat)}

    ch_taxid_bam = ch_taxid_bam.mix(SAMTOOLS_SORT_PASS.out.bam)
        .join (ch_mapped_reads, by: [0])
        .filter { _meta, _bam, _mapped, pass -> pass }
        .map { meta, bam, _mapped, _pass -> [meta, bam] }
    ch_taxid_bai = ch_taxid_bai.mix(SAMTOOLS_INDEX_PASS.out.index)
        .join(ch_mapped_reads, by:[0])
        .filter { _meta, _bai, _mapped, pass -> pass }
        .map { meta, bai, _mapped, _pass -> [meta, bai] }

    // Prepare individual FASTA files for each taxID with the number of mapped reads less than params.min_read_counts
    ch_blast_input = ch_accession_taxid_with_meta.fail
        .join( input_bam, by: 0 ) // Join by meta (index 0)
        .map { meta, accession_list, taxid, organism, bam, bam_index ->
            // Create new meta with taxid and organism information
            def new_meta = meta + [taxid: taxid, organism: organism, accessions: accession_list.flatten()]
            return [ new_meta, bam, bam_index]
        }

    // FASTA files will be used as BLAST input, bam file will be used in IGV
    SAMTOOLS_VIEW_FAIL(ch_blast_input, [[],[],[]], [[],[]], [[],[]], [] )

    ch_fail_subset = SAMTOOLS_VIEW_FAIL.out.bam
        .map { meta, bam -> [meta.subMap(meta.keySet() - 'accessions'), bam ] }

    SAMTOOLS_SORT_FAIL(ch_fail_subset, [[],[],[]], 'bai')
    SAMTOOLS_INDEX_FAIL(SAMTOOLS_SORT_FAIL.out.bam)

    SAMTOOLS_FASTA(SAMTOOLS_SORT_FAIL.out.bam, false)

    // Combine pathogen-specific BAM/BAI files from both pass and fail branches
    ch_pathogen_bam_bai = channel.empty()
    ch_pathogen_bam_bai = ch_pathogen_bam_bai.mix(ch_taxid_bam)
        .join(ch_taxid_bai, by:0)
        .mix(SAMTOOLS_SORT_FAIL.out.bam
                .join(SAMTOOLS_INDEX_FAIL.out.index, by:0)
        )

    // Generate pathogen-specific coverage and depth
    SAMTOOLS_COVERAGE (ch_pathogen_bam_bai, [[],[],[]])
    SAMTOOLS_DEPTH (ch_pathogen_bam_bai, [[],[]])

    ch_coverage_plot_input = channel.empty()
    ch_coverage_plot_input = ch_coverage_plot_input
        .mix(SAMTOOLS_DEPTH.out.tsv.filter {_meta, depth_file -> depth_file.size() > 0 })
        .join(SAMTOOLS_COVERAGE.out.coverage, by:0)

    COVERAGE_PLOT(ch_coverage_plot_input)


    emit:
    taxid_bam       = ch_taxid_bam
    taxid_bai       = ch_taxid_bai
    taxid_fasta     = SAMTOOLS_FASTA.out.fasta
    taxid_bam_fail  = SAMTOOLS_SORT_FAIL.out.bam
    taxid_bai_fail  = SAMTOOLS_INDEX_FAIL.out.index
    coverage        = SAMTOOLS_COVERAGE.out.coverage
    depth           = SAMTOOLS_DEPTH.out.tsv
    coverage_plot   = COVERAGE_PLOT.out.png
}
