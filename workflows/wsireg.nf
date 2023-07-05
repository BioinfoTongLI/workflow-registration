#!/usr/bin/env/ nextflow

// Copyright (C) 2023 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tifs_in = [:]
params.out_dir = "./out/"
params.enhance = 0

params.enable_conda = false

params.generate_fake_anchor = true
params.max_n_worker = 30
params.ref_ch = "DAPI" // or dapi
params.ref_cycle = 0

DEBUG=true

process wsi_reg {
    debug DEBUG

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioinfotongli/registration:wsireg':
        'bioinfotongli/registration:wsireg'}"
    /*containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"*/
    /*publishDir params.out_dir + "/wsireg_output", mode: 'copy'*/
    storeDir params.out_dir + "/wsireg_output"

    input:
    tuple val(meta), path(ref), val(meta_moving), path(moving)
    val(target_ch_index)

    output:
    tuple val(meta_moving), path("${wsireg_output}.tiff"), emit: registered_tif
    tuple val(meta_moving), path("${wsireg_output}_transformations.json"), emit: tmats, optional: true

    script:
    meta['stem'] = ref.baseName
    ref_ind = meta["stem"]
    moving_ind = meta_moving["index"]
    wsireg_output = "${ref_ind}/${moving_ind}-${moving_ind}_to_${ref_ind}"
    """
    mkdir $ref_ind # somehow wsireg won't create this folder on its own
    wsireg_wrapper_2_images.py \
        -ref_ind ${ref_ind} \
        -ref_img ${ref} \
        -moving_ind ${moving_ind} \
        -moving_img ${moving} \
        -target_ch_index ${target_ch_index} \
    """
}
