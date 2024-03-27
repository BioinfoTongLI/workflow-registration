VERSION = '0.5.0'
DOCKERHUB_IMAGE = "bioinfotongli/raw2ometiff:${VERSION}"

process BIOINFOTONGLI_RAW2OMETIFF {
    tag "$meta.id"
    label 'process_medium'

    /*conda (params.enable_conda ? "-c ome raw2ometiff=${VERSION}" : null)*/
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${DOCKERHUB_IMAGE}":
        "${DOCKERHUB_IMAGE}" }"
    publishDir params.out_dir , mode: 'copy'

    input:
    tuple val(meta), path(ome_zarr)

    output:
    tuple val(meta), path("${meta.id}.ome.tif"), emit: ome_tif
    path "raw2ometiff_versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (meta['id'] == null) {
        meta['id'] = ome_zarr.baseName
    }
    def args = task.ext.args ?: ''
    """
    raw2ometiff \\
        --max_workers=$task.cpus \\
        $args \\
        $ome_zarr \\
        ${meta.id}.ome.tif

    cat <<-END_VERSIONS > raw2ometiff_versions.yml
    "${task.process}":
        bioinfotongli: \$(echo \$(raw2ometiff --version 2>&1) | sed 's/^.*raw2ometiff //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
