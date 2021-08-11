#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tifs_in = ""
params.out_dir = ""

params.generate_fake_anchor = true
params.double_feature_reg = false
params.tilesize = 1000
params.max_n_worker = 12
params.ref_ch = "dapi"

include { Feature_based_registration; fake_anchor_chs; Second_register} from workflow.projectDir + '/registration.nf'
include { bf2raw; raw2bf} from workflow.projectDir + '/image-convert/main.nf'

workflow {
    Channel.fromPath(params.ome_tifs_in)
        .map{it: file(it)}
        .collect()
        .set{ome_tif_paths}
    Feature_based_registration(ome_tif_paths, params.ref_ch)
    /*fake_anchor_ch(feature_based_registration.out)*/
    Feature_based_registration.out.view()
    Second_register(Feature_based_registration.out, params.ref_ch)
    bf2raw(Second_register.out)
    raw2bf(bf2raw.out)
}
