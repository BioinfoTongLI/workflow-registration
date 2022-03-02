#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tifs_in = ""
params.out_dir = ""

params.tilesize = 1000
params.max_n_worker = 12
params.ref_ch = "DAPI"
params.ref_cycle = 0

params.singularity_img = "/nfs/cellgeni/singularity/images/registration-v0.0.1.sif"


process Feature_based_registration {
    tag "${images}"
    container params.singularity_img
    containerOptions "--containall"
    /*publishDir params.out_dir, mode:"copy"*/
    /*storeDir params.out_dir + "/first_reg"*/
    /*clusterOptions = '-pe smp 10 -l virtual_free=64G,h_rt=30:00:00'*/

    cpus params.max_n_worker
    memory "50G"

    input:
    path(images)
    val ref_ch

    output:
    path("out.tif"), emit: feature_reg_tif

    script:
    """
    python /feature_reg/reg.py -i ${images} -o ./ -r "${params.ref_cycle}" -c "${ref_ch}" -n ${params.max_n_worker} --tile_size ${params.tilesize}
    """
}


process OpticalFlow_register {
    tag "${tif}"
    echo true

    container params.singularity_img
    containerOptions "--containall"
    /*storeDir params.out_dir + "/second_reg"*/

    cpus params.max_n_worker
    memory "150G"

    input:
    path(tif)
    val ref_ch

    output:
    path('*opt_flow_registered.tif')

    script:
    """
    python /opt_flow_reg/opt_flow_reg.py -i "${tif}" -c "${ref_ch}" -o ./ -n ${params.max_n_worker} --tile_size ${params.tilesize} --overlap 100  --method rlof
    """
}


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


workflow {
    Channel.fromPath(params.ome_tifs_in)
        .map{it: file(it)}
        .collect()
        .set{ome_tif_paths}
    Feature_based_registration(ome_tif_paths, params.ref_ch)
    OpticalFlow_register(Feature_based_registration.out, params.ref_ch)
    bf2raw(OpticalFlow_register.out)
    raw2bf(bf2raw.out)
}
