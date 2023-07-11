
process BIOINFOTONGLI_MICROALIGNER {
    tag "$meta.id"
    label 'process_large'

    /*conda "YOUR-TOOL-HERE"*/ // Do not support this yet
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioinfotongli/microaligner:develop':
        'bioinfotongli/microaligner:develop' }"

    input:
    tuple val(meta), path(config), path(images)

    output:
    tuple val(meta), path(f"${prefix}.ome.tif"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    microaligner \\
        ${config}
        $args \\
        -o ${prefix}.ome.tif \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(microaligner --version 2>&1) | sed 's/^.*microaligner //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(microaligner --version 2>&1) | sed 's/^.*microaligner //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
