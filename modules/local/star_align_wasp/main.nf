// Local fork of nf-core/modules star/align @ 6d46786420b4d7bc88eba026eb389c0c5535d120.
//
// Why it's forked: STAR's WASP filter (--waspOutputMode SAMtag) needs a
// per-sample VCF of heterozygous variants staged into the task directory,
// and the upstream module has no input slot for one. WASP re-maps each read
// with its alleles swapped and discards any read whose placement changes —
// that's what removes reference-allele mapping bias, which otherwise shows
// up downstream as fake allele-specific expression in phASER.
//
// The ONLY differences from upstream are:
//   - renamed STAR_ALIGN -> STAR_ALIGN_WASP
//   - an added `tuple val(meta4), path(vcf)` input
//   - the $wasp line in the script block
// Keep it that way: when upstream changes, re-copy the module and re-apply
// those three edits rather than letting this drift into a rewrite.
// The unmodified upstream module is still installed at
// modules/nf-core/star/align/ so `nf-core modules update` can show the diff.

process STAR_ALIGN_WASP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data' :
        'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4' }"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(vcf) // heterozygous variants for WASP; pass [[],[]] to align without WASP
    val star_ignore_sjdbgtf

    output:
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed "s/STAR_//"'), emit: versions_star, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('gawk'), eval("gawk --version | sed -n '1s/GNU Awk \\([0-9.]*\\).*/\\1/p'"), emit: versions_gawk, topic: versions

    tuple val(meta), path('*d.out.bam')                              , optional:true, emit: bam
    tuple val(meta), path("${prefix}.sortedByCoord.out.bam")         , optional:true, emit: bam_sorted
    tuple val(meta), path("${prefix}.Aligned.sortedByCoord.out.bam") , optional:true, emit: bam_sorted_aligned
    tuple val(meta), path('*toTranscriptome.out.bam')                , optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam')                 , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')                               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                                   , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')                            , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')                  , optional:true, emit: read_per_gene_tab
    tuple val(meta), path('*.out.junction')                          , optional:true, emit: junction
    tuple val(meta), path('*.out.sam')                               , optional:true, emit: sam
    tuple val(meta), path('*.wig')                                   , optional:true, emit: wig
    tuple val(meta), path('*.bg')                                    , optional:true, emit: bedgraph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads1 = []
    def reads2 = []
    meta.single_end ? [reads].flatten().each{ read -> reads1 << read} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }
    def ignore_gtf      = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' 'SM:$prefix'"
    def out_sam_type    = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    // STAR reads the VCF with plain file I/O — it must be uncompressed, and
    // must carry genotypes (GT) so STAR can tell which sites are het.
    // Het genotypes are all it needs: phased (0|1) and unphased (0/1) VCFs
    // give identical results here.
    // The vA/vG/vW SAM attributes that carry the WASP verdict per read are
    // set through ext.args (see conf/modules.config), not here.
    def wasp            = vcf ? "--waspOutputMode SAMtag --varVCFfile ${vcf}" : ''
    mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $ignore_gtf \\
        $attrRG \\
        $wasp \\
        $args

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}Xd.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg
    """
}
