process EXTRACTDIAMONDREADS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebc86c3786bdc70fe31e9fac77ab4d05af3d27664a259298a831d1994fcb022d/data'
        : 'community.wave.seqera.io/library/seqkit_pigz_python:3ff8d721e22b2875'}"

    input:
    val taxid
    val evalue_threshold
    tuple val(meta), path(tsv), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz"), optional: true, emit: extracted_diamond_reads
    tuple val("${task.process}"), val("python"), eval("python --version | sed -e 's/Python //g'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val("pigz"), eval("pigz --version 2>&1 | sed -e 's/pigz //g'"), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single_end_flag = meta.single_end ? '--single_end' : ''

    """
    extractdiamondreads.py \\
        --taxid ${taxid} \\
        --tsv ${tsv} \\
        --evalue ${evalue_threshold} \\
        --prefix ${prefix} \\
        ${single_end_flag} \\
        --fastq ${fastq.join(' ')} # joins the list of file paths into a space-separated string

    # Compress the resulting fastq files
    pigz -p ${task.cpus} *.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzip_reads_command = meta.single_end ?
        "gzip ${prefix}.extracted_diamond_read.fastq" :
        "gzip ${prefix}.extracted_diamond_read_1.fastq; gzip ${prefix}.extracted_diamond_read_2.fastq"
    """
    if [ "${meta.single_end}" == 'true' ]; then
        echo '@stub_read' > ${prefix}.extracted_diamond_read.fastq
        echo 'ACGT' >> ${prefix}.extracted_diamond_read.fastq
        echo '+' >> ${prefix}.extracted_diamond_read.fastq
        echo '!!!!' >> ${prefix}.extracted_diamond_read.fastq
        $gzip_reads_command
    else
        echo '@stub_read/1' > ${prefix}.extracted_diamond_read_1.fastq
        echo 'ACGT' >> ${prefix}.extracted_diamond_read_1.fastq
        echo '+' >> ${prefix}.extracted_diamond_read_1.fastq
        echo '!!!!' >> ${prefix}.extracted_diamond_read_1.fastq
        echo '@stub_read/2' > ${prefix}.extracted_diamond_read_2.fastq
        echo 'TGCA' >> ${prefix}.extracted_diamond_read_2.fastq
        echo '+' >> ${prefix}.extracted_diamond_read_2.fastq
        echo '!!!!' >> ${prefix}.extracted_diamond_read_2.fastq
        $gzip_reads_command
    fi
    """
}
