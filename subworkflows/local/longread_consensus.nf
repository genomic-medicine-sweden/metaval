//
// Consensus calling for long reads
//

include { SAMTOOLS_CONSENSUS } from '../../modules/nf-core/samtools/consensus/main'
include { MEDAKA } from '../../modules/nf-core/medaka/main'
include { GUNZIP } from '../../modules/nf-core/gunzip/main'

workflow LONGREAD_CONSENSUS {
    take:
    bam         // channel: [ val(meta), path(bam) ]
    reference   // channel: [ val(meta), path(reference) ]

    main:
    ch_versions = Channel.empty()

    if ( params.longread_consensus_tool == 'medaka' ) {
        input_medaka  = bam.combine(reference).map{ meta_bam, bam, meta_ref, reference -> [ meta_bam, bam, reference ]}
        MEDAKA ( input_medaka )
        ch_consensus_gz = MEDAKA.out.assembly
        ch_versions = ch_versions.mix(MEDAKA.out.versions)

        GUNZIP ( ch_consensus_gz )
        ch_consensus = GUNZIP.out.gunzip
        ch_versions = ch_versions.mix(GUNZIP.out.versions)
    } else if ( params.longread_consensus_tool == 'samtools' ) {
        // Combine bam and reference channels
        SAMTOOLS_CONSENSUS ( bam )
        ch_consensus = SAMTOOLS_CONSENSUS.out.fasta
        ch_versions = ch_versions.mix(SAMTOOLS_CONSENSUS.out.versions)
    }

    emit:
    consensus = ch_consensus // channel: [ val(meta), path(consensus) ]
    versions = ch_versions   // channel: [ versions.yml ]
}
