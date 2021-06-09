#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tifs_in = ""
params.out_dir = ""

params.generate_fake_anchor = true
params.double_feature_reg = false
params.tilesize = 1000
params.max_n_worker = 12

include { Feature_based_registration; fake_anchor_chs; Second_register} from workflow.projectDir + '/registration.nf'

workflow {
    Channel.fromPath(params.ome_tifs_in)
        .map{it: file(it)}
        .collect()
        .set{ome_tif_paths}
    Feature_based_registration(ome_tif_paths)
    /*fake_anchor_ch(feature_based_registration.out)*/
    Feature_based_registration.out.view()
    Second_register(Feature_based_registration.out)
}
