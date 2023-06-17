#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

/*params.ome_tifs_in = "*/
params.feature_reg_yaml = ""
params.optflow_reg_yaml = ""
params.rna_molecule_radius = 3
params.out_dir = "./out/"
params.peak_ch_index = 3
params.target_ch_index = 0
params.trackpy_percentile = 70
params.search_range = [5]
params.tracking = false
params.enhance = 0

params.enable_conda = false

params.generate_fake_anchor = true
params.double_feature_reg = false
params.tilesize = 1000
params.ref_ch = "DAPI" // or dapi
params.stem = "20220511_hindlimb"
params.sif_folder = null 
params.ref_cycle = 0


//include { BIOINFOTONGLI_BIOFORMATS2RAW as bf2raw} from '/lustre/scratch126/cellgen/team283/tl10/modules/modules/bioinfotongli/bioformats2raw/main.nf' addParams(
//    enable_conda:false,
//    publish:false,
//    store:true,
//    out_dir:params.out_dir
//)

process Feature_based_registration {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "microaligner.sif":
        'microaligner:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /scratch,/mnt':'-v /scratch'}"

    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    path(config_file)
    /*path(images)*/
    /*val ref_ch*/
    /*val ref_cycle*/

    output:
    path("*feature_reg_result_stack.tif"), emit: feature_reg_tif
    /*path("tmats.tsv"), emit: feature_reg_tmat*/

    script:
    """
    microaligner ${config_file}
    """
}
/*python /opt/feature_reg/reg.py -i ${images} -o ./ -r ${ref_cycle} -c "${ref_ch}" -n ${params.max_n_worker} --tile_size ${params.tilesize}*/
/*featurereg.py -imgs $images*/


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

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "microaligner.sif":
        'microaligner:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /scratch,/mnt':'-v /scratch'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    path(tif)
    path(config_file_for_optflow)

    output:
    tuple val(meta), path("*_optflow_reg_result_stack.tif")

    script:
    meta = ["stem":params.stem]
    """
    microaligner ${config_file_for_optflow}
    """
}
/*python /opt/opt_flow_reg/opt_flow_reg.py -i "${tif}" -c "${ref_ch}" -o ./ -n ${params.max_n_worker} --tile_size ${params.tilesize} --overlap 100 # --method rlof*/
/*mv out_opt_flow_registered.tif ${stem}_opt_flow_registered.tif*/

process Wsireg {
    debug false

    label "default"
    label "large_mem"
    /*label "medium_mem"*/

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "wsireg.sif":
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    /*containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"*/

    storeDir params.out_dir + "/wsireg_output"

    input:
    /*tuple path(ref), path(moving)*/
    tuple val(ref_ind), path(ref_image), val(moving_ind), path(moving_img)
    val(target_ch_index)

    output:
    tuple val(moving_ind), path("${ref_ind}/${moving_ind}-${moving_ind}_to_${ref_ind}_registered.ome.tiff"), emit: registered_tif
    tuple val(moving_ind), path("${ref_ind}/${moving_ind}-${moving_ind}_to_${ref_ind}_transformations.json"), emit: tmats, optional: true

    script:
    """
    mkdir $ref_ind # somehow wsireg won't create this folder on its own
    wsireg_wrapper_2_images.py \
        -ref_ind "${ref_ind}" \
        -ref_img ${ref_image} \
        -moving_ind "${moving_ind}" \
        -moving_img ${moving_img} \
        -target_ch_index ${target_ch_index} \
    """
}

process Call_peaks {

    label "default"
    label "large_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "wsireg.sif":
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    storeDir params.out_dir + "/peaks_after_registration"

    input:
    tuple val(ind), path(one_cycle_image)
    val(spot_radius)
    val(percentile)
    val(peak_ch_index)

    output:
    path("${stem}_peaks_radius_${spot_radius}_percentile_${percentile}_with_profile.tsv")

    script:
    stem = file(one_cycle_image).baseName
    """
    call_fake_peaks.py -stem $stem -spot_radius ${spot_radius} -frame ${ind} -percentile ${percentile} -peak_ch_index ${peak_ch_index} $one_cycle_image
    """
}


process Whitehat_tiff {

    /*label "cellgeni_a100"*/
    /*label "gpu_normal"*/
    label "imaging"
    label "large_mem"
    /*label "medium_mem"*/

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "wsireg.sif":
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    storeDir params.out_dir + "/whitehated_tiffs"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"

    cpus 10

    input:
    tuple val(ind), path(img)
    val(spot_radius)
    val(enhance)

    output:
    tuple val(ind), path("${stem}_whitehat_enhanced_radius_${spot_radius}.ome.tif")

    script:
    stem = img.baseName
    """
    whitehat_tiff.py -stem "$stem" -tiff ${img} -spot_radius ${spot_radius} -enhance ${enhance}
    """
}

process Track_peaks {

    label "small_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "wsireg.sif":
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    storeDir params.out_dir + "/tracked_peaks"

    input:
    path(tsvs)
    val(search_range)

    output:
    path("tracked_peaks_search_range_${search_range}.tsv")

    script:
    """
    track_molecules.py -search_range ${search_range} ${tsvs}
    """
}


process Filter_tracks {

    label "small_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "wsireg.sif":
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    storeDir params.out_dir + "/filtered_tracks"

    input:
    path(tracks)

    output:
    path("${stem}_filtered.tsv")

    script:
    stem = tracks.baseName
    """
    filter_tracks.py -stem ${stem} -tracks ${tracks}
    """
}

process Extract_intensity_profile_of_anchor {

    label "large_mem"
    queue "imaging"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "wsireg.sif":
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    storeDir params.out_dir + "/anchor_peak_profile"

    input:
    tuple path(anchor_peaks), val(frame), path(tiff)

    output:
    path("anchor_profile_of_frame_${frame}.tsv")

    script:
    """
    extract_intensity_profile_of_anchor.py -anchor_peaks ${anchor_peaks} -frame ${frame} \
        -tiff ${tiff}
    """
}

process Prepare_profile_for_decoding {
    debug true

    /*label "medium_mem"*/
    /*label "huge_mem"*/
    label "large_mem"
    cpus 1

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_folder + "wsireg.sif":
        'gitlab-registry.internal.sanger.ac.uk/tl10/workflow-registration:wsireg'}"
    storeDir params.out_dir + "/profiles_for_decoding"

    input:
    path(peak_profiles)

    output:
    path("anchor_profile_for_decoding.npy")

    script:
    """
    prepare_profile_for_decoding.py ${peak_profiles}
    """
}

/*ome_tif_paths = Channel.from(params.ome_tifs_in)*/

/*paired_ome_tif_paths = Channel.fromPath(params.ome_tifs_in[0]).combine(Channel.fromPath(params.ome_tifs_in[1, 2, 3, 4, 5, 6]))*/
/*paired_ome_tif_paths.view()*/

workflow {
    Feature_based_registration(channel.fromPath(params.feature_reg_yaml))
    /*Feature_based_registration(ome_tif_paths.map{it: file(it[1])}.collect(), params.ref_ch, params.ref_cycle)*/
    /*fake_anchor_ch(feature_based_registration.out)*/
    OpticalFlow_register(Feature_based_registration.out, channel.fromPath(params.optflow_reg_yaml))
    //bf2raw(OpticalFlow_register.out)
}

workflow Featurereg {
    take:
        images
    main:
        Feature_based_registration(images,
            channel.from(params.ref_ch),
            channel.from(params.ref_cycle)
        )
    emit:
        Feature_based_registration.out
}

workflow Align_peaks {
    Whitehat_tiff(ome_tif_paths, params.rna_molecule_radius, params.enhance)

    all_cycles = Whitehat_tiff.out.branch{
        ref_cycle: it[0] == 0
        coding_cycles: it
    }
    cycles_to_register = all_cycles.ref_cycle.combine(all_cycles.coding_cycles)

    Wsireg(cycles_to_register, params.target_ch_index)

    if (params.tracking) {
        Call_peaks(Wsireg.out.registered_tif, params.rna_molecule_radius, params.trackpy_percentile, params.peak_ch_index)
        /*Call_peaks.out.collect().view()*/
        Track_peaks(Call_peaks.out.collect(), channel.from(params.search_range))

        Filter_tracks(Track_peaks.out)
    } else {
        Call_peaks(all_cycles.ref_cycle, params.rna_molecule_radius, params.trackpy_percentile, params.peak_ch_index)
        /*Call_peaks.out.view()*/
        Extract_intensity_profile_of_anchor(Call_peaks.out.combine(Wsireg.out.registered_tif))
        Prepare_profile_for_decoding(Extract_intensity_profile_of_anchor.out.collect())
    }
}
