process MEDAKA_PARALLEL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/13ffeff27cc2e3ba985a2d436f9352ff6ea17c8366d81a5c00a50b3a253aa9a5/data'
        : 'community.wave.seqera.io/library/medaka_seqkit:0cb6c6d35f0ad431'}"

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*_sorted.fasta"), emit: assembly
    tuple val("${task.process}"), val("medaka"), eval("medaka --version 2>&1 | sed 's/medaka //g'"), emit: versions_medaka, topic: versions
    tuple val("${task.process}"), val("seqkit"), eval("seqkit | sed '3!d; s/Version: //'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Export prefix to bash so it always exists
    prefix='${prefix}'

    mkdir -p tmp
    export TMPDIR=./tmp

    # Safe decompression of input assembly
    assembly_path="${assembly}"
    if [[ "\$assembly_path" == *.gz ]]; then
        # Decompress to a temporary file without touching original
        gunzip -c "\$assembly_path" > "\${prefix}_assembly.fa"
        assembly="\${prefix}_assembly.fa"
    fi

    mini_align \\
        -i ${reads} \\
        -r \$assembly \\
        -m \\
        -t ${task.cpus} \\
        -p ${prefix}_calls_to_draft \\
        ${args1}

    # In medaka >= 2.0 this step is medaka inference, in earlier versions it is consensus
    mkdir inference
    # Start with the largest contigs, they probably take longest
            # Medaka can do with 2 threads and may need some extra for IO

    sort -nrk2 \${assembly}.fai \\
        | cut -f1 | xargs -P ${task.cpus} \\
        -n1 \\
        -I{} \\
        medaka inference ${prefix}_calls_to_draft.bam \\
            inference/{}.hdf \\
            --region {} \\
            --threads 2 \\
            ${args2}

    # In medaka >= 2.0 this step is medaka sequence, in earlier versions it is stitch
    medaka sequence \\
        --threads ${task.cpus} \\
        ${args3} \\
        inference/*.hdf \$assembly ${prefix}.fa

    # Sort the consensus by reads ID
    seqkit sort -n ${prefix}.fa > ${prefix}_sorted.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sorted.fasta
    """
}
