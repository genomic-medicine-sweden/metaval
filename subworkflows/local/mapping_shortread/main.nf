//
// Mapping
//

include { BOWTIE2_BUILD                             } from '../../../modules/nf-core/bowtie2/build'
include { FASTQ_ALIGN_BOWTIE2                       } from '../../nf-core/fastq_align_bowtie2'
include { SAMTOOLS_FAIDX                            } from '../../../modules/nf-core/samtools/faidx'
include { PIGZ_UNCOMPRESS                           } from '../../../modules/nf-core/pigz/uncompress'
include { SAMTOOLS_SORT                             } from '../../../modules/nf-core/samtools/sort'
include { getFlagstatMappedReads                    } from '../utils_nfcore_metaval_pipeline'
include { SAMTOOLS_COVERAGE                         } from '../../../modules/nf-core/samtools/coverage'
include { SAMTOOLS_DEPTH                            } from '../../../modules/nf-core/samtools/depth'
include { COVERAGE_PLOT                             } from '../../../modules/local/coverage_plot'

workflow MAPPING_SHORTREAD {
    take:
    ch_reads_reference // [ [ meta ], [ reads ], [ reference ] ]
    generate_coverage // Boolean: whether to generate coverage and depth files for the BAM files

    main:
    ch_multiqc_files = channel.empty()
    ch_coverage      = channel.empty()
    ch_depth         = channel.empty()
    ch_coverage_plot = channel.empty()

    // Build the bowtie2 index
    BOWTIE2_BUILD (
        ch_reads_reference.map { meta, _reads, ref -> [ meta, ref ] }
    )
    ch_bowtie2_index = BOWTIE2_BUILD.out.index

    // Build the reference index fai
    PIGZ_UNCOMPRESS (
        ch_reads_reference.map { meta, _reads, ref -> [ meta, ref ] }
    )
    ch_ref_uncompressed = PIGZ_UNCOMPRESS.out.file
        .map { meta, ref -> [meta, ref, []]}

    SAMTOOLS_FAIDX ( ch_ref_uncompressed, false )

    // Join the uncompressed reference and reference index
    ch_ref_fai = ch_ref_uncompressed
        .map {meta, ref, _empty -> [meta, ref]}
        .join(SAMTOOLS_FAIDX.out.fai, by:0)

    // Join the reads, minimap2 index, reference and reference index
    ch_reads_with_index = ch_reads_reference
        .map { meta, reads, _ref -> [ meta, reads ] }
        .join(ch_bowtie2_index, by: 0)
        .join(ch_ref_fai, by:0)
        .multiMap { meta, reads, bowtie2_index, ref, fai ->
            reads: [meta, reads]
            index: [meta, bowtie2_index]
            fasta_fai: [meta, ref, fai]
        }

    // Mapping
    FASTQ_ALIGN_BOWTIE2 (
        ch_reads_with_index.reads,
        ch_reads_with_index.index,
        false,                          // save unaligned
        false,                          // sort bam
        ch_reads_with_index.fasta_fai
    )

    // Remove empty bam files
    ch_flagstat = channel.empty()
    ch_flagstat = ch_flagstat.mix(FASTQ_ALIGN_BOWTIE2.out.flagstat)
        .map { meta, flagstat -> [meta] + getFlagstatMappedReads(flagstat)}

    ch_bam_mapped = channel.empty()
    ch_bam_mapped = ch_bam_mapped.mix(FASTQ_ALIGN_BOWTIE2.out.bam)
        .join (ch_flagstat, by: [0])
        .map { meta, bam, _mapped, pass -> if (pass) [meta, bam] }

    SAMTOOLS_SORT(ch_bam_mapped, [[],[],[]], 'bai')

    // Generate samtools stats for coverage and depth
    if (generate_coverage) {
        ch_bam_bai = channel.empty()
        ch_bam_bai = ch_bam_bai.mix(SAMTOOLS_SORT.out.bam)
            .join(SAMTOOLS_SORT.out.index, by:0)

        SAMTOOLS_COVERAGE (ch_bam_bai, [[],[],[]])
        SAMTOOLS_DEPTH (ch_bam_bai, [[],[]])

        ch_coverage_plot_input = channel.empty()
        ch_coverage_plot_input = ch_coverage_plot_input.mix(SAMTOOLS_DEPTH.out.tsv.filter {_meta, depth_file -> depth_file.size() > 0 })
            .join(SAMTOOLS_COVERAGE.out.coverage, by:0)

        COVERAGE_PLOT (ch_coverage_plot_input)

        ch_coverage = SAMTOOLS_COVERAGE.out.coverage
        ch_depth = SAMTOOLS_DEPTH.out.tsv
        ch_coverage_plot = COVERAGE_PLOT.out.png

    }


    ch_multiqc_files = ch_multiqc_files.mix ( FASTQ_ALIGN_BOWTIE2.out.flagstat.collect{ _meta, flagstat_file -> flagstat_file }.ifEmpty([]) )

    emit:
    index         = BOWTIE2_BUILD.out.index           // channel: [ val(meta), [ index ] ]
    bam           = SAMTOOLS_SORT.out.bam             // channel: [ val(meta), [ bam ] ]
    bai           = SAMTOOLS_SORT.out.index           // channel: [ val(meta), [ bai ] ]
    coverage      = ch_coverage                       // channel: [ val(meta), [ coverage ] ]
    depth         = ch_depth                          // channel: [ val(meta), [ depth ] ]
    coverage_plot = ch_coverage_plot                  // channel: [ val(meta), [ coverage_plot ] ]
    mqc           = ch_multiqc_files
}
