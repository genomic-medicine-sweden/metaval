process FLAG_TAXPASTA {
    tag "${meta1.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'biocontainers/pandas:2.2.1'}"

    input:
    tuple val(meta1), path(sample_taxpasta, stageAs: 'sample_taxpasta.tsv'), val(meta2), path(ntc_taxpasta, stageAs: 'ntc_taxpasta.tsv')

    output:
    tuple val(meta1), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val("python"), eval("python --version | sed -e 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val("pandas"), eval("python -c \"import pandas; print(pandas.__version__)\""), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta1.id}"
    def ntc_cmd = ntc_taxpasta ? "--ntc-taxpasta ${ntc_taxpasta} --ntc-name ${meta2.id}" : ''

    """
    flag_taxpasta.py \\
        --sample-taxpasta ${sample_taxpasta} \\
        --sample-name ${meta1.id} \\
        ${ntc_cmd} \\
        --output ${prefix}.tsv \\
        ${args}
    """
}
