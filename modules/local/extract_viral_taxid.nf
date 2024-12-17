process EXTRACT_VIRAL_TAXID {

    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.8.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.2--h9ee0642_1':
        'biocontainers/seqkit:2.8.2--h9ee0642_1' }"

    input:
    tuple val(meta), path(taxpasta_standardised_profile)
    tuple val(meta), path(report) // classification report

    output:
    tuple val(meta), path("*viral_taxids.tsv"), optional:true, emit: viral_taxid
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.tool}"

    """
    if grep -qi "virus" $taxpasta_standardised_profile; then
        grep -i "virus" $taxpasta_standardised_profile | cut -f 1 > taxpasta_viral_taxid.txt
        if [[ "${meta.tool}" == "kraken2" || "${meta.tool}" == "centrifuge" ]]; then
            awk -F'\t' '\$3 != 0 {print \$5}' ${report} > detected_taxid.txt
            grep -F -w -f taxpasta_viral_taxid.txt detected_taxid.txt > ${prefix}_viral_taxids.tsv
        elif [[ "${meta.tool}" == "diamond" ]]; then
            cut -f 2 ${report} | uniq > detected_taxid.txt
            grep -F -w -f taxpasta_viral_taxid.txt detected_taxid.txt | uniq > ${prefix}_viral_taxids.tsv
        fi
    else
        echo "No viral taxids found." > "no_viral_taxid.txt"
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
