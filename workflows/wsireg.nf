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
params.sif_folder = "/lustre/scratch126/cellgen/team283/imaging_sifs/"
params.ref_cycle = 0

DEBUG=true

process wsireg {
    debug DEBUG

    label "default"
    label "large_mem"
    /*label "medium_mem"*/

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "wsireg.sif":
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    /*containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"*/

    storeDir params.out_dir + "/wsireg_output"

    input:
    /*tuple path(ref), path(moving)*/
    tuple val(ref_ind), path(ref_image), val(moving_ind), path(moving_img)
    val(target_ch_index)

    output:
    tuple val(moving_ind), path("${ref_ind}/${moving_ind}-${moving_ind}_to_${ref_ind}_registered.ome.tiff"), emit: registered_tif
    tuple val(moving_ind), path("${ref_ind}/${moving_ind}-${moving_ind}_to_${ref_ind}_transformations.json"), emit: tmats, optional: true

    script:
    """
    mkdir $ref_ind # somehow wsireg won't create this folder on its own
    wsireg_wrapper_2_images.py \
        -ref_ind "${ref_ind}" \
        -ref_img ${ref_image} \
        -moving_ind "${moving_ind}" \
        -moving_img ${moving_img} \
        -target_ch_index ${target_ch_index} \
    """
}
