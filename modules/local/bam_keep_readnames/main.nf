// Local module — no nf-core equivalent exists.
//
// Keeps only the records whose read name survived UMI deduplication.
//
// Why this exists: UMI-tools dedup can only run on the GENOME BAM. It
// decides two reads are PCR duplicates by comparing their UMI *and* their
// alignment coordinate, and a transcriptome BAM's coordinates are per
// transcript — the same fragment appears once per compatible isoform, at a
// different position each time. Running dedup there independently would both
// double-count and disagree with the genome BAM about which reads are real.
//
// So dedup happens ONCE, on the genome BAM, and this module propagates that
// single verdict: pull the surviving read names out of the deduplicated
// genome BAM, then keep exactly those reads in the transcriptome BAM. Every
// transcript record of a surviving fragment is kept (Salmon needs all of
// them to assign the read across isoforms); every record of a fragment
// dedup threw away is dropped.
//
// `samtools view -N` needs samtools >= 1.12; the container below is 1.24.

process BAM_KEEP_READNAMES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam)       // BAM to filter (the transcriptome BAM)
    tuple val(meta2), path(name_bam) // BAM whose read names define what to keep (the deduplicated genome BAM)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("*.readnames.txt"), optional: true, emit: readnames
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.filtered"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    # -T . keeps sort's spill files inside the task dir; a whole run's worth
    # of read names is easily a few GB and would otherwise fill /tmp.
    samtools view -@ ${task.cpus} ${name_bam} \\
        | cut -f1 \\
        | LC_ALL=C sort -u -T . \\
        > ${prefix}.readnames.txt

    samtools view \\
        -@ ${task.cpus} \\
        -b \\
        -N ${prefix}.readnames.txt \\
        -o ${prefix}.bam \\
        ${args} \\
        ${bam}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.filtered"
    """
    touch ${prefix}.bam
    touch ${prefix}.readnames.txt
    """
}
