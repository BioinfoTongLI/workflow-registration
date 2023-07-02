#!/usr/bin/env/ nextflow

// Copyright (C) 2023 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.enable_conda = false

params.feature_reg_yaml = ""
params.optflow_reg_yaml = ""
params.max_n_worker = 30
params.out_dir = "./out/"

params.sif_folder = "/lustre/scratch126/cellgen/team283/imaging_sifs/"

DEBUG=true

include { BIOINFOTONGLI_TO_OME_TIFF } from '../subworkflows/sanger/bioinfotongli/to_ome_tiff/main' addParams(
    enable_conda:false,
    publish:false,
    store:true,
    out_dir:params.out_dir
)

process Feature_based_registration {
    debug DEBUG

    label 'default'
    label 'large_mem'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "microaligner.sif":
        'microaligner:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    cpus params.max_n_worker

    input:
    path(config_file)

    output:
    path("*feature_reg_result_stack.tif"), emit: feature_reg_tif
    /*path("tmats.tsv"), emit: feature_reg_tmat*/

    script:
    """
    microaligner ${config_file}
    """
}

process OpticalFlow_register {
    debug DEBUG

    label 'default'
    label 'large_mem'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "microaligner.sif":
        'microaligner:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    path(tif)
    path(config_file_for_optflow)

    output:
    tuple val(meta), path("*_optflow_reg_result_stack.tif")

    script:
    meta = ["stem":params.stem]
    """
    microaligner ${config_file_for_optflow}
    """
}


workflow micro_aligner {
    Feature_based_registration(channel.fromPath(params.feature_reg_yaml))
    /*fake_anchor_ch(feature_based_registration.out)*/
    OpticalFlow_register(Feature_based_registration.out, channel.fromPath(params.optflow_reg_yaml))
    BIOINFOTONGLI_TO_OME_TIFF(OpticalFlow_register.out)
}
