#!/usr/bin/env/ nextflow

// Copyright (C) 2023 Tong LI <tongli.bioinfo@proton.me>

nextflow.enable.dsl=2

include { micro_aligner } from './workflows/microaligner'
include { itk_reg } from './workflows/itk'

params.images = [
    [['index': 0], '/nfs/team283_imaging/SM_BRA/playground_Tong/Susanna_Jimmy_hiplex/20230608_registration/SM_0037_l/SM_0037_l-Cycle_1.ome.tif'],
    [['index': 1], '/nfs/team283_imaging/SM_BRA/playground_Tong/Susanna_Jimmy_hiplex/20230608_registration/SM_0037_l/SM_0037_l-Cycle_2.ome.tif'],
    [['index': 2], '/nfs/team283_imaging/SM_BRA/playground_Tong/Susanna_Jimmy_hiplex/20230608_registration/SM_0037_l/SM_0037_l-Cycle_3.ome.tif'],
]
params.ref_cycle = 1

channel.from(params.images)
        .branch{ it ->
            ref: it[0]['index'] == params.ref_cycle
            movings: it[0]['index'] != params.ref_cycle
        }
        .set{images}
images.ref.combine(images.movings).view()


workflow run_micro_aligner {
    micro_aligner()
}

workflow run_itk_reg {
    itk_reg(
        images.ref.combine(images.movings)
    )
}
