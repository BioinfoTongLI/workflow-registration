#!/usr/bin/env/ nextflow

// Copyright (C) 2023 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.out_dir = "./out/"
params.images = "/nfs/team283_imaging/image0001.tif"

DEBUG=true

process QCAlignment {
    debug DEBUG

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioinfotongli/registration:qc':
        'bioinfotongli/registration:qc'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    publishDir params.out_dir, mode : 'copy'

    input:
    tuple val(meta), path(tif)

    output:
    path(out_folder)

    script:
    out_folder = "${meta.stem}_QC"
    """
    AlignmentQC.py \
        --filepath ${tif} \
        --output_folder ${out_folder} \
        --Ncyc ${meta.Ncyc} \
        --DapiCh ${meta.DapiCh}
    """
}
