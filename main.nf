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
params.stem = "20220511_hindlimb"
params.sif_folder = "/lustre/scratch117/cellgen/team283/tl10/sifs/"

/*
 * bf2raw: The bioformats2raw application converts the input image file to
 * an intermediate directory of tiles in the output directory.
 */
process bf2raw {
    debug true
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
    debug true
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
    debug true
    /*container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"*/
    container params.sif_folder + "feature_reg.sif"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    path(images)
    val ref_ch

    output:
    path("out.tif"), emit: feature_reg_tif

    script:
    """
    python /opt/feature_reg/reg.py -i ${images} -o ./ -r 0 -c "${ref_ch}" -n ${params.max_n_worker} --tile_size ${params.tilesize}
    """
}


process fake_anchor_chs {
    debug true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/generate_fake_anchors:latest"
    containerOptions "-B ${baseDir}:/code:ro"
    /*storeDir params.out_dir + "/first_reg"*/
    /*publishDir params.out_dir, mode:"copy"*/

    input:
    file(ome_tif)

    output:
    path("*anchors.ome.tif"), emit: tif_with_anchor

    script:
    if (params.generate_fake_anchor){
        """
        python /code/generate_fake_anchors.py -ome_tif ${ome_tif} -known_anchor "c01 Alexa 647"
        """
    } else {
        """
        #rename it
        """
    }
}

process OpticalFlow_register {
    debug true
    /*container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"*/
    container params.sif_folder + "opt_flow_reg.sif"
    storeDir params.out_dir

    input:
    path tif
    val ref_ch
    val stem

    output:
    path("${stem}_opt_flow_registered.tif")

    script:
    """
    python /opt/opt_flow_reg/opt_flow_reg.py -i "${tif}" -c "${ref_ch}" -o ./ -n ${params.max_n_worker} --tile_size ${params.tilesize} --overlap 100 # --method rlof
    mv out_opt_flow_registered.tif ${stem}_opt_flow_registered.tif
    """
}

process wsireg {
    debug true

    /*label "default"*/

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'my.sif':
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"

    storeDir params.out_dir + "/wsireg_output"

    input:
    path(images)

    output:
    /*tuple val(stem), path("${stem}.out")*/
    path("*.ome.tif")

    script:
    /*stem = file(images).baseName*/
    stem = "Test"
    """
    mkdir $stem
    wsireg_wrapper.py -stem $stem ${images}
    """
}

workflow {
    ome_tif_paths = Channel.fromPath(params.ome_tifs_in)
        .map{it: file(it)}
        .collect()
    Feature_based_registration(ome_tif_paths, params.ref_ch)
    /*fake_anchor_ch(feature_based_registration.out)*/
    OpticalFlow_register(Feature_based_registration.out, params.ref_ch, params.stem)
    bf2raw(OpticalFlow_register.out)
    raw2bf(bf2raw.out)
}

workflow Wsireg {
    ome_tif_paths = Channel.fromPath(params.ome_tifs_in)
        .map{it: file(it)}
        .collect()
    wsireg(ome_tif_paths)
}
