process METAVAL_REPORT {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b45d221e51a26945c244afa9bd126a6757289be51b963721ce9f25f7c8662c38/data'
        : 'community.wave.seqera.io/library/biopython_jinja2_pandas_python:bf9cf8457c0990de'}"

    input:
    val ticket
    path samplesheet
    path flagged_taxpasta , stageAs: 'flagged_taxpasta/*'
    path reads            , stageAs: 'reads/*'
    path blastn           , stageAs: 'blastn/*'
    path blastx           , stageAs: 'blastx/*'
    path coverage_tables  , stageAs: 'coverage_tables/*'
    path coverage_plots   , stageAs: 'coverage_plots/*'

    output:
    path("*.html"), emit: metaval_report
    tuple val("${task.process}"), val("python"), eval("python --version | sed -e 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val("jinja2"), eval("python -c \"import jinja2; print(jinja2.__version__)\""), emit: versions_jinja2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'metaval_report'
    def pipeline_version = workflow.manifest.version ?: 'dev'

    """
    static_report.py \\
        --samplesheet ${samplesheet} \\
        --flagged-dir ./flagged_taxpasta \\
        --reads-dir ./reads \\
        --blastn-dir ./blastn \\
        --blastx-dir ./blastx \\
        --coverage-dir ./coverage_tables \\
        --coverage-plots-dir ./coverage_plots \\
        --version ${pipeline_version} \\
        --ticket ${ticket} \\
        --output ${prefix}.html \\
        ${args}
    """
}
