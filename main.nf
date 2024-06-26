#!/usr/bin/env/ nextflow

// Copyright (C) 2023 Tong LI <tongli.bioinfo@proton.me>

include { micro_aligner } from './workflows/microaligner'
include { wsi_reg } from './workflows/wsireg'
include { QCAlignment } from './workflows/alignQC'

params.images = [
    [['index': 0], 'image1'],
    [['index': 1], 'image2'],
]
params.ref_cycle = 0
params.target_ch_index = 0

channel.from(params.images)
        .branch{ it ->
            ref: it[0]['index'] == params.ref_cycle
            movings: it[0]['index'] != params.ref_cycle
        }
        .set{images}
images.ref.combine(images.movings)
    .set{registeration_pairs}


process stack {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioinfotongli/registration:wsireg':
        'bioinfotongli/registration:wsireg'}"
    /*containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"*/
    publishDir params.out_dir, mode:"copy"

    input:
    tuple val(meta), path(ref)
    path(registered_tifs)

    output:
    tuple val(meta), path(out), emit: hyperstack

    script:
    meta['stem'] = ref.baseName
    out = "${meta['stem']}_registered_all_cycle_stack.ome.tif"
    """
    stack_tifs.py \
        -ref ${ref} \
        -ref_index ${meta['index']} \
        -out ${out} \
        ${registered_tifs}
    """
}


workflow run_micro_aligner {
    micro_aligner()
    // QCAlignment(micro_aligner.out)
}

workflow run_wsireg {
    wsi_reg(registeration_pairs, params.target_ch_index)
    wsi_reg.out.registered_tif
        .map{ it -> file(it[1]) }
        .collect(sort:true)
        .set{registered_tifs}
    stack(images.ref, registered_tifs)
    QCAlignment(stack.out.hyperstack)
}
