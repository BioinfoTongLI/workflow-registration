#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tifs_in = ""
params.out_dir = ""

params.double_feature_reg = false
params.tilesize = 1000
params.max_n_worker = 12
params.ref_ch = "DAPI" // or dapi


/*
 * bf2raw: The bioformats2raw application converts the input image file to
 * an intermediate directory of tiles in the output directory.
 */
process bf2raw {
    echo true
    conda "-c ome bioformats2raw"
    /*storeDir params.out_dir + "/raws"*/
    /*publishDir params.out_dir, mode:"copy"*/

    input:
    path(img)

    output:
    tuple val(stem), file("${stem}")

    script:
    stem = img.baseName
    """
    bioformats2raw --max_workers ${params.max_n_worker} --resolutions 7 --no-hcs $img "${stem}"
    """
}


/*
 * raw2bf: The raw2ometiff application converts a directory of tiles to
 * an OME-TIFF pyramid.
 */
process raw2bf {
    echo true
    conda "-c ome raw2ometiff"
    /*storeDir params.out_dir + "/ome_with_pyramids"*/
    publishDir params.out_dir, mode:"copy"

    input:
    tuple val(stem), path(zarr)

    output:
    path("${stem}.tif")

    script:
    """
    raw2ometiff --max_workers ${params.max_n_worker} ${zarr} "${stem}.tif"
    """
}


process Feature_based_registration {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"
    /*containerOptions "--cpus=${params.max_n_worker}"*/
    /*publishDir params.out_dir, mode:"copy"*/
    /*storeDir params.out_dir + "/first_reg"*/

    cpus params.max_n_worker

    input:
    val(ref_ch)
    path(images)

    output:
    path("out.tif"), emit: feature_reg_tif

    script:
    """
    python /feature_reg/reg.py -i ${images} -o ./ -r 0 -c "${ref_ch}" -n ${params.max_n_worker} --tile_size ${params.tilesize}
    """
}

process optflow_register {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"
    containerOptions "--cpus=${params.max_n_worker}"
    /*storeDir params.out_dir + "/second_reg"*/

    input:
    tuple val(ref_ch), path(tif)

    output:
    path('*opt_flow_registered.tif')

    script:
    """
    python /opt_flow_reg/opt_flow_reg.py -i "${tif}" -c "${ref_ch}" -o ./ -n ${params.max_n_worker} --tile_size ${params.tilesize} --overlap 100  --method rlof
    """
}


process Local_Register {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"
    /*containerOptions "--cpus=${params.max_n_worker}"*/
    storeDir params.out_dir + "/local_registered"

    input:
    tuple val(ref_ch), path(zarr)

    output:
    path('*local_registered.tif')

    script:
    """
    python ${projectDir}/local_reg.py --zarr_root "${zarr}"
    #-c "${ref_ch}" -o ./ -n ${params.max_n_worker} --tile_size ${params.tilesize} --overlap 100  --method rlof
    """
}


workflow {
    Channel.fromPath(params.ome_tifs_in)
        .map{it: file(it)}
        .collect()
        .set{ome_tif_paths}
    Feature_based_registration(params.ref_ch, ome_tif_paths)
    /*Feature_based_registration.out.view()*/
    bf2raw(Feature_based_registration.out.feature_reg_tif)
    Local_Register(bf2raw.out)
    /*Local_Register()*/
    /*optflow_register(Feature_based_registration.out, params.ref_ch)*/
    /*bf2raw(Second_register.out)*/
    /*raw2bf(bf2raw.out)*/
}
