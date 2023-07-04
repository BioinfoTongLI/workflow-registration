#!/usr/bin/env/ nextflow

// Copyright (C) 2023 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2


params.out_dir = "./out/"

params.sif_folder = "/lustre/scratch126/cellgen/team283/imaging_sifs/"

DEBUG=true

include { BIOINFOTONGLI_TO_OME_TIFF } from '../subworkflows/sanger/bioinfotongli/to_ome_tiff/main' addParams(
    publish:false,
    store:true,
    out_dir:params.out_dir
)

process itk_reg {
    debug DEBUG

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioinfotongli/registration:itk':
        'bioinfotongli/registration:itk'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    publishDir params.out_dir, mode:"copy"

    input:
    tuple val(meta), path(ref), val(meta_moving), path(moving)

    output:
    tuple val(meta), path("*itk_reg_stack.tif"), emit: itk_reg_tif
    path("tmats.tsv"), emit: itk_reg_tmat, optional:true

    script:
    """
    itk_2_images.py -ref ${ref} -moving ${moving}
    """
}

/*workflow itk_reg {*/
    /*itk_reg(channel.fromPath(params.images))*/
    /*BIOINFOTONGLI_TO_OME_TIFF(itk_reg.out)*/
/*}*/
