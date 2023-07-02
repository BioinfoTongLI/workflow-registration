//
// Convert any microscopy image to ome.tiff
//

/*params.align_options          = [:]*/
params.bf2raw_options = [:]

include { BIOINFOTONGLI_RAW2OMETIFF as raw2ometiff } from '../../../../modules/sanger/bioinfotongli/raw2ometiff/main' addParams( options: params.bf2raw_options    )
include { BIOINFOTONGLI_BIOFORMATS2RAW as bf2raw } from '../../../../modules/sanger/bioinfotongli/bioformats2raw/main' addParams( sort_options: params.samtools_sort_options, )

workflow BIOINFOTONGLI_TO_OME_TIFF {
    take:
    images // channel: [ val(meta), file(image) ]

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with Bowtie2
    //
    bf2raw ( reads, index )
    ch_versions = ch_versions.mix(bf2raw.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    raw2ometiff ( bf2raw.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    emit:
    bam_orig          = bf2raw.out.bam          // channel: [ val(meta), bam   ]
    log_out           = bf2raw.out.log          // channel: [ val(meta), log   ]

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]

    versions       = ch_versions                      // channel: [ versions.yml ]
}
