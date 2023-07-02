VERSION = '0.5.0'

process BIOINFOTONGLI_BIOFORMATS2RAW {
    tag "${meta.id}"
    label 'process_medium'

    cpus task.cpus

    conda (params.enable_conda ? "-c ome bioformats2raw==${VERSION}" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "openmicroscopy/bioformats2raw:${VERSION}":
        "openmicroscopy/bioformats2raw:${VERSION}" }"

    input:
    tuple val(meta), path(img)

    output:
    path "{$stem}.zarr", emit: zarr
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def stem = meta['id'] ?: img.baseName
    def args = task.ext.args ?: ''

    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    """
    /opt/conda/bin/bioformats2raw \\
        $args \\
        --max_workers $task.cpus \\
        $img \\
        "${stem}.zarr"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioinfotongli: \$(echo \$(bioformats2raw --version 2>&1) | sed 's/^.*bioformats2raw //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
