#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tifs_in = ""
params.out_dir = ""

params.generate_fake_anchor = true
params.double_feature_reg = false
params.tilesize = 1000
params.max_n_worker = 12
params.ref_ch = "DAPI" // or dapi

include { Feature_based_registration; fake_anchor_chs; Second_register} from workflow.projectDir + '/registration.nf'

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
    /*fake_anchor_ch(feature_based_registration.out)*/
    /*Feature_based_registration.out.view()*/
    Second_register(Feature_based_registration.out, params.ref_ch)
    bf2raw(Second_register.out)
    raw2bf(bf2raw.out)
}
