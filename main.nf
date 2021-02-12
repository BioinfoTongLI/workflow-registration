#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

params.ome_tifs_in = [""]
params.out_dir = "/nfs/team283_imaging/JSP_HSS/playground_Tong/human_brain_158_20x_reg_double_feature_reg_small_tile_4_4_DAPI_4_4_fake_anchor/"

params.generate_fake_anchor = true
params.double_feature_reg = true
params.max_n_worker = 12

Channel.fromPath(params.ome_tifs_in)
    .map{it: file(it)}
    .collect()
    .set{ome_tif_paths}

process feature_based_registration {
    /*container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"*/
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:small_tiles"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/first_reg"

    input:
    file images from ome_tif_paths

    output:
    file "out.tif" into feature_based_regs

    script:
    """
    python /image_registrator/reg.py -i ${images} -o ./ -r 0 -c "DAPI" -n ${params.max_n_worker}
    """
}

process fake_anchor_chs {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/generate_fake_anchors:latest"
    containerOptions "-v ${baseDir}:/code:ro"
    storeDir params.out_dir + "/first_reg"
    /*publishDir params.out_dir, mode:"copy"*/

    input:
    file ome_tif from feature_based_regs

    output:
    file "*anchors.ome.tif" into tif_with_anchor, tif_with_anchor_for_pyramid

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

process Second_register {
    echo true
    /*container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:grid_10"*/
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:small_tiles"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/second_reg"

    input:
    file tif from tif_with_anchor

    output:
    file expected_out into opt_registered

    script:
    if (params.double_feature_reg){
        expected_out = 'out.tif'
        """
        python /image_registrator/reg.py -i ${tif} -o ./ -r 0 -c "anchor" -n ${params.max_n_worker} --stack
        """
    } else {
        expected_out = '*registered.ome.tif'
        """
        python /opt_flow_reg/opt_flow_reg.py -i "${tif}" -c "anchor" -o ./ -n ${params.max_n_worker}
        """
    }
}

process convert_to_n5 {
    /*echo true*/
    container "gitlab-registry.internal.sanger.ac.uk/tl10/img-bftools"
    /*storeDir params.out_dir + "/raws"*/
    /*publishDir params.out_dir + "/tmp", mode:"copy"*/

    input:
    file img from opt_registered

    output:
    tuple val(stem), file("${stem}") into raws

    script:
    stem = img.baseName
    """
    /bf2raw/bioformats2raw-0.2.6/bin/bioformats2raw --dimension-order XYZCT --max_workers ${params.max_n_worker} --resolutions 7 --tile_width 512 --tile_height 512 $img "${stem}"
    """
}


process n5_to_ome_tiff {
    /*echo true*/
    container "gitlab-registry.internal.sanger.ac.uk/tl10/img-bftools"
    /*storeDir params.out_dir + "/second_reg_with_pyramid"*/
    publishDir params.out_dir + "/second_reg_with_pyramid", mode:"copy"

    input:
    tuple val(stem), file(zarr) from raws

    output:
    tuple val(stem), path("${stem}.ome.tif") into ome_tiffs

    script:
    """
    /raw2tif/raw2ometiff-0.2.8/bin/raw2ometiff --max_workers ${params.max_n_worker} --debug "${zarr}" "${stem}.ome.tif"
    """
}
