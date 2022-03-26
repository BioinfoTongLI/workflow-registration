#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.tilesize = 2000

process Feature_based_registration {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"
    /*containerOptions "--cpus=${params.max_n_worker}"*/
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    path(images)
    val ref_ch

    output:
    path("out.tif"), emit: feature_reg_tif

    script:
    """
    python /feature_reg/reg.py -i ${images} -o ./ -r 0 -c "${ref_ch}" -n ${params.max_n_worker} --tile_size ${params.tilesize}
    """
}


process fake_anchor_chs {
    echo true
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

process Second_register {
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"
    /*containerOptions "--cpus=${params.max_n_worker}"*/
    storeDir params.out_dir

    input:
    path(tif)
    val ref_ch

    output:
    path(expected_out)

    script:
    if (params.double_feature_reg){
        expected_out = 'out.tif'
        """
        python /feature_reg/reg.py -i ${tif} -o ./ -r 0 -c "anchor" -n ${params.max_n_worker} --stack --tile_size ${params.tilesize}
        """
    } else {
        expected_out = '*opt_flow_registered.tif'
        """
        python /opt_flow_reg/opt_flow_reg.py -i "${tif}" -c "${ref_ch}" -o ./ -n ${params.max_n_worker} --tile_size ${params.tilesize} --overlap 100  --method rlof
        """
    }
}
