//
// Prepare an indivudual BAM/FASTA file for each pathogen with mapped reads
//

include { SUBSET_BAM                                        } from '../../modules/local/subset_bam'
include { SAMTOOLS_SORT                                     } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                                    } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS                                 } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_FILTER     } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FASTA                                    } from '../../modules/nf-core/samtools/fasta/main'
include { GUNZIP as GUNZIP_SINGLE                           } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PAIRED                           } from '../../modules/nf-core/gunzip/main'

workflow TAXID_BAM_FASTA {
    take:
    bam
    bai
    accession2taxid
    min_read_counts

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    input_bam =  bam.combine( bai,by: 0 )
    SAMTOOLS_IDXSTATS( input_bam )
    ch_accession = SAMTOOLS_IDXSTATS.out.idxstats
        .map { it[1] }
        .splitCsv( header: false,sep:"\t" )
        .filter { it -> it[0]!= "*" }

    ch_versions.mix( SAMTOOLS_IDXSTATS.out.versions.first() )

    // Load accession2taxid.map
    ch_accession2taxidmap = accession2taxid.splitCsv( header: false,sep:"\t" )

    ch_accession_taxid = ch_accession2taxidmap
        .join( ch_accession )
        .filter { it -> it[3] != "0" }
        .map { [ it[0], it[1] ] }
        .groupTuple( by: 1 )

    ch_samtools_view = ch_accession_taxid
        .combine(input_bam)
        .map {accession_list, taxid, meta, bam, bam_index ->
            def new_meta = meta.clone()
            new_meta.taxid = taxid
            return [ new_meta, bam, bam_index, accession_list ]
        }
        .multiMap {
            meta, bam, bam_index, accession_list ->
                bam:  [meta, bam, bam_index]
                accession: accession_list.flatten()
        }

    // Individual BAM file for each taxID
    SUBSET_BAM ( ch_samtools_view.bam, ch_samtools_view.accession )
    ch_versions  = ch_versions.mix( SUBSET_BAM.out.versions.first() )

    SAMTOOLS_SORT ( SUBSET_BAM.out.bam, [[],[]] )
    ch_versions  = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions  = ch_versions.mix( SAMTOOLS_INDEX.out.versions.first() )


    // Count reads in each BAM file
    SAMTOOLS_IDXSTATS_FILTER ( SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai) )
    ch_versions = ch_versions.mix( SAMTOOLS_IDXSTATS_FILTER.out.versions.first() )
    // Process IDXSTATS output to get total read count
    ch_bam_counts = SAMTOOLS_IDXSTATS_FILTER.out.idxstats
        .map { meta, stats_file ->
            def total_reads = 0
            stats_file.eachLine { line ->
                def fields = line.split('\t')
                total_reads += fields[2].toInteger() + fields[3].toInteger()
            }
            [ meta, total_reads ]
        }
    // Filter BAM files based on read count: only call
    ch_filtered_bam = SAMTOOLS_SORT.out.bam
        .join(ch_bam_counts)
        .filter { meta, bam, count -> count >= min_read_counts }
        .map { meta, bam, count -> [meta, bam] }

    // Convert filtered BAM files to FASTA format
    SAMTOOLS_FASTA ( SAMTOOLS_SORT.out.bam, false )
    ch_versions = ch_versions.mix( SAMTOOLS_FASTA.out.versions.first() )

    // Separate single-end and paired-end FASTA files
    ch_fasta_single = SAMTOOLS_FASTA.out.fasta.filter { meta, fasta -> meta.single_end }
    ch_fasta_paired = SAMTOOLS_FASTA.out.fasta.filter { meta, fasta -> !meta.single_end }

    // Unzip single-end fasta files
    GUNZIP_SINGLE ( ch_fasta_single )
    ch_versions  = ch_versions.mix( GUNZIP_SINGLE.out.versions )

    // Unzip paired-end fasta files
    GUNZIP_PAIRED ( ch_fasta_paired.transpose() )
    ch_versions  = ch_versions.mix( GUNZIP_PAIRED.out.versions )

    ch_unzipped_fasta = GUNZIP_SINGLE.out.gunzip.mix(GUNZIP_PAIRED.out.gunzip.groupTuple())

    emit:
    versions        = ch_versions
    taxid_bam       = ch_filtered_bam
    taxid_fasta     = ch_unzipped_fasta
    taxid_accession = ch_samtools_view.accession

}
