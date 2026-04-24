process COVERAGE_PLOT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/matplotlib:3.5.1'
        : 'biocontainers/matplotlib:3.5.1'}"

    input:
    tuple val(meta), path(depth_file), path(coverage_file)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val("${task.process}"), val("python"), eval("python --version | sed -e 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    coverage_plot.py \\
        --depth_file ${depth_file} \\
        --coverage_file ${coverage_file} \\
        --output ${prefix}.png \\
        ${args}
    """
}
