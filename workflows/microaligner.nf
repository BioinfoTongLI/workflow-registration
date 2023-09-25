#!/usr/bin/env/ nextflow

// Copyright (C) 2023 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.enable_conda = false

params.feature_reg_yaml = ""
params.optflow_reg_yaml = ""
params.out_dir = "./out/"

DEBUG=true

include { BIOINFOTONGLI_MICROALIGNER as featreg; BIOINFOTONGLI_MICROALIGNER as optreg } from '../modules/sanger/bioinfotongli/microaligner/main'

include { BIOINFOTONGLI_TO_OME_TIFF } from '../subworkflows/sanger/bioinfotongli/to_ome_tiff/main' addParams(
    enable_conda:false,
    publish:false,
    store:true,
    out_dir:params.out_dir
)

VERSION = 'v1.0.6'


process Feature_based_registration {
    debug DEBUG

    label 'default'
    label 'large_mem'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "bioinfotongli/microaligner:${VERSION}":
        "bioinfotongli/microaligner:${VERSION}"}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    path(config_file)

    output:
    path("*feature_reg_result_stack.tif"), emit: feature_reg_tif
    /*path("tmats.tsv"), emit: feature_reg_tmat*/

    script:
    """
    microaligner --config ${config_file}
    """
}

process OpticalFlow_register {
    debug DEBUG

    label 'default'
    label 'large_mem'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "bioinfotongli/microaligner:${VERSION}":
        "bioinfotongli/microaligner:${VERSION}"}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    path(tif)
    path(config_file_for_optflow)

    output:
    tuple val(meta), path("*_optflow_reg_result_stack.tif")

    script:
    meta = ["stem": tif.stem, "ext": tif.ext]
    """
    microaligner --config ${config_file_for_optflow}
    """
}


workflow micro_aligner {
    main:
        Feature_based_registration(channel.fromPath(params.feature_reg_yaml))
        /*fake_anchor_ch(feature_based_registration.out)*/
        OpticalFlow_register(Feature_based_registration.out, channel.fromPath(params.optflow_reg_yaml))
        BIOINFOTONGLI_TO_OME_TIFF(OpticalFlow_register.out)
    emit:
        BIOINFOTONGLI_TO_OME_TIFF.out.ome_tif
}

workflow {
    main:
        featreg(channel.fromPath(params.images))
        optreg(featfeatreg.out.registered_image)
    emit:
        optreg.out.ome_tif
}
