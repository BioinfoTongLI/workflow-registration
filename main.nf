#!/usr/bin/env/ nextflow

// Copyright (C) 2020 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

include { micro_aligner } from './workflows/microaligner'

workflow run_micro_aligner {
    micro_aligner()
}
