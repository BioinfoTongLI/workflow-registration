#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

params.ome_tifs_in = ["/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-02T09_58_03-Measurement 8b/JSP_HSS_OB10037_W9-BRA_Nucleus_RCPs_Meas8b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-02T18_23_21-Measurement 13b/JSP_HSS_OB10037_W9-BRA_Nucleus_b1A_b1C_b1T_b1G_Meas13b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-05T14_22_34-Measurement 20b/JSP_HSS_OB10037_W9-BRA_Nucleus_b2A_b2C_b2T_b2G_Meas20b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/201006_122322-V__2020-10-06T12_34_07-Measurement 1b/JSP_HSS_OB10037_W9-BRA_Nucleus_b3A_b3C_b3T_b3G_Meas1b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-07T10_08_00-Measurement 35b/JSP_HSS_OB10037_W9-BRA_Nucleus_b4A_b4C_b4T_b4G_Meas35b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-08T10_41_11-Measurement 44b/JSP_HSS_OB10037_W9-BRA_Nucleus_b5A_b5C_b5T_b5G_Meas44b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-10T00_01_32-Measurement 55b/JSP_HSS_OB10037_W9-BRA_Nucleus_b6A_b6C_b6T_b6G_Meas55b_A3_F1T0.ome.tif"]
/*params.out_dir = "/nfs/team283_imaging/JSP_HSS/playground_Tong/human_brain_158_opt_flow_registered_20x_non_normalization/"*/
params.out_dir = "/nfs/team283_imaging/JSP_HSS/playground_Tong/human_brain_158_opt_flow_registered_20x_normalized/"
/*params.out_dir = "/nfs/team283_imaging/JSP_HSS/playground_Tong/human_brain_158_20x_reg_double_feature_reg/"*/
params.out_dir = "/nfs/team283_imaging/JSP_HSS/playground_Tong/human_brain_158_20x_reg_double_feature_reg_small_tile/"

/*params.ome_tifs_in = ["/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-02T18_23_21-Measurement 13b/JSP_HSS_OB10037_W9-BRA_Nucleus_b1A_b1C_b1T_b1G_Meas13b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-05T14_22_34-Measurement 20b/JSP_HSS_OB10037_W9-BRA_Nucleus_b2A_b2C_b2T_b2G_Meas20b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/201006_122322-V__2020-10-06T12_34_07-Measurement 1b/JSP_HSS_OB10037_W9-BRA_Nucleus_b3A_b3C_b3T_b3G_Meas1b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-07T10_08_00-Measurement 35b/JSP_HSS_OB10037_W9-BRA_Nucleus_b4A_b4C_b4T_b4G_Meas35b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-08T10_41_11-Measurement 44b/JSP_HSS_OB10037_W9-BRA_Nucleus_b5A_b5C_b5T_b5G_Meas44b_A3_F1T0.ome.tif", "/nfs/team283_imaging/0HarmonyStitched/JSP_HSS/JSP_HSS_OB10037__2020-10-10T00_01_32-Measurement 55b/JSP_HSS_OB10037_W9-BRA_Nucleus_b6A_b6C_b6T_b6G_Meas55b_A3_F1T0.ome.tif"]*/
/*params.out_dir = "/nfs/team283_imaging/JSP_HSS/playground_Tong/human_brain_158_20x_cyc2_7_reg/"*/

params.generate_fake_anchor = true
params.double_feature_reg = true

Channel.fromPath(params.ome_tifs_in)
    .map{it: file(it)}
    .collect()
    .set{ome_tif_paths}

process feature_based_registration {
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/first_reg"

    input:
    file images from ome_tif_paths

    output:
    file "out.tif" into feature_based_regs

    script:
    """
    python /image_registrator/reg.py -i ${images} -o ./ -r 0 -c "DAPI" -n 12
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
    container "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:latest"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/second_reg"

    input:
    file tif from tif_with_anchor

    output:
    file "*.tif" into opt_registered

    script:
    if (params.double_feature_reg){
        """
        python /image_registrator/reg.py -i ${tif} -o ./ -r 0 -c "anchor" -n 12 --stack
        """
    } else {
        """
        python /opt_flow_reg/opt_flow_reg.py -i "${tif}" -c "anchor" -o ./ -n 14
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
    /bf2raw/bioformats2raw-0.2.6/bin/bioformats2raw --dimension-order XYZCT --max_workers 15 --resolutions 7 --tile_width 512 --tile_height 512 $img "${stem}"
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
    export _JAVA_OPTIONS="-Xmx128g"
    /raw2tif/raw2ometiff-0.2.8/bin/raw2ometiff --max_workers 12 --debug "${zarr}" "${stem}.ome.tif"
    """
}
