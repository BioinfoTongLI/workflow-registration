#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tifs_in = "
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
params.max_n_worker = 30
params.ref_ch = "DAPI" // or dapi
params.sif_folder = "/lustre/scratch126/cellgen/team283/imaging_sifs/"
params.ref_cycle = 0


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

workflow align_peaks {
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
