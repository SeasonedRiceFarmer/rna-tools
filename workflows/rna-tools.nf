/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { BBMAP_BBSPLIT          } from '../modules/nf-core/bbmap/bbsplit/main'
include { KRAKEN2_KRAKEN2        } from '../modules/nf-core/kraken2/kraken2/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SORTMERNA                } from '../modules/nf-core/sortmerna/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rna-tools_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNA_TOOLS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // Absurdly overcommented because I'm a high school chud who
    // lowk dunno anything about this lmao.. for learning

    // ── STEP 1: FASTQC ──────────────────────────────────────
    // Reports raw read quality (per-base scores, GC content,
    // adapter contamination, duplication) — diagnostic only,
    // does not modify reads.

    FASTQC(ch_samplesheet)

    // Collect FastQC's zip reports for the final MultiQC summary.
    // MultiQC reads zips directly, so no metadata needed here —
    // just the file paths. Note: mix combines emissions 
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    // Build fastp's input tuple: [meta, reads, adapter_fasta].
    // Empty list = no adapter file, fastp auto-detects adapters.
    // TODO: expose adapter_fasta as an optional pipeline param
    ch_fastp_input = ch_samplesheet.map { meta, reads -> [ meta, reads, [] ] }

    // TODO: Optional view for learning, marked for removal.
    ch_fastp_input.view()

    // ── STEP 2: FASTP ──────────────────────────────────────
    // Independent to FastQC, this is the trimming step for raw reads.
    // What it does:
    //   - Clips adapter contamination off read ends; purely
    //     statistical if no adapter_fasta is given.
    //   - trims low-quality bases, typically from the 3' end where
    //     quality tends to degrade first
    //   - drops whole reads that are too short or too low-quality
    //     after trimming to be usable
    // Outputs: cleaned FASTQ files (fed to STAR downstream) +
    // its own HTML/JSON report (pre vs post-trim stats), which
    // also gets collected into ch_multiqc_files.
    // TODO: Add an option for "trim galore!" or other tools. 

    FASTP(
        ch_fastp_input,
        false, // discard_trimmed_pass: keep reads that failed trimming, don't silently drop them
        false, // save_trimmed_fail: don't save a separate file of failed reads
        false, // save_merged: don't merge overlapping PE reads into one combined read
    )

    // What you get after trimming, fastp emitted to reads. Tuple of reads.
    ch_trimmed_reads = FASTP.out.reads

     // TODO: Optional view for learning, marked for removal.
    ch_trimmed_reads.view()

    // Updates the channel so that the channel can later be given to MultiQC. (Fastp JSON report)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{ _meta, file -> file })

    // ── STEP 3: FASTQC ──────────────────────────────────────
    // FastQC again to verify Fastp successfully trimmed
    // for higher quality.

    // Tag the meta.id with a "_trimmed" suffix so the output
    // files (and thus filenames) are distinguishable from the
    // raw FastQC run above — MultiQC groups by filename pattern,
    // and both runs would otherwise collide in one section.
    ch_fastqc_trimmed_input = ch_trimmed_reads.map { meta, reads -> [ meta + [id: "${meta.id}_trimmed"], reads ] }

    FASTQC_TRIMMED(ch_fastqc_trimmed_input)

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.map{ _meta, file -> file })

    // ── STEP 4: Kraken2 (optional) ──────────────────────────
    // Classifies trimmed reads against a broad reference DB to see
    // which non-human taxa (bacteria/fungi/virus) are actually
    // present in each sample. Its report doesn't remove anything by
    // itself — it just tells STEP 5 (BBSplit) which extra genomes
    // are worth screening against, on top of whatever's already in
    // --bbsplit_fasta_list. Off by default: needs a large pre-built
    // Kraken2 DB (e.g. PlusPF).

    // Per-sample taxids that clear both thresholds (rel. abundance
    // of total reads + absolute read count); empty channel when
    // Kraken2 didn't run, which .collect() below still resolves to [].
    ch_kraken2_taxids = channel.empty()

    if (params.run_kraken2) {
        if (!params.kraken2_db) {
            error "Please provide --kraken2_db (path to a Kraken2 database) when --run_kraken2 is true."
        }

        KRAKEN2_KRAKEN2(
            ch_trimmed_reads,
            file(params.kraken2_db, checkIfExists: true),
            false, // save_output_fastqs: we only need the report, not classified/unclassified fastqs
            false, // save_reads_assignment: don't need the per-read taxid assignment file
        )

        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.map{ _meta, file -> file })

        // TODO: Optional view for learning, marked for removal.
        KRAKEN2_KRAKEN2.out.report.view { meta, report -> "[kraken2] report for ${meta.id}: ${report}" }

        // Kraken2's report.txt is tab-separated with no header:
        // percent_of_reads, reads_in_clade, reads_direct, rank_code, taxid, name
        // The report has one row per taxonomic RANK, not just species — a
        // human read rolls up through root/Eukaryota/Metazoa/.../Homo sapiens,
        // and every one of those ancestor rows trivially clears any abundance
        // threshold when most reads are human. Only rank_code 'S'-prefixed rows
        // (species, or S1/S2 strain-level) are actual organisms with a genome
        // worth screening for — and Homo sapiens (9606) is the host, not a
        // contaminant, so it's excluded explicitly.
        ch_kraken2_taxids = KRAKEN2_KRAKEN2.out.report
            .flatMap { _meta, report ->
                report.readLines()
                    .findAll { line ->
                        def cols = line.tokenize('\t')
                        def percent     = cols[0].trim() as Double
                        def reads_clade = cols[1].trim() as Integer
                        def rank_code   = cols[3].trim()
                        def taxid       = cols[4].trim()
                        rank_code.startsWith('S') && taxid != '9606' &&
                            percent >= (params.kraken2_min_rel_abundance * 100) && reads_clade >= params.kraken2_min_reads
                    }
                    .collect { line -> line.tokenize('\t')[4].trim() }
            }

        // TODO: Optional view for learning, marked for removal.
        ch_kraken2_taxids.view { taxid -> "[kraken2] taxid cleared both thresholds: ${taxid}" }
    }

    // ── STEP 5: BBSplit (optional) ──────────────────────────
    // Screens trimmed reads against extra "contaminant" genomes
    // (bacteria/fungi/virus) alongside the primary human genome,
    // and keeps only the reads that best match the human genome.
    // Contaminant genomes come from two places, merged together:
    //   - --bbsplit_fasta_list: a static, user-curated CSV
    //   - Kraken2 (if enabled): taxa actually detected in this run,
    //     resolved to reference fastas via --kraken2_taxid_lookup
    // One shared index is built from the union of both across the
    // whole run and reused for every sample — index building is the
    // expensive part, and this keeps every sample's stats comparable
    // in MultiQC. Off by default: --perform_bbsplit false.

    ch_filtered_reads = ch_trimmed_reads

    if (params.perform_bbsplit) {
        if (!params.bbsplit_fasta_list && !params.run_kraken2) {
            error "Please provide --bbsplit_fasta_list and/or enable --run_kraken2 (contaminant genomes for BBSplit to screen against) when --perform_bbsplit is true."
        }
        if (!params.fasta) {
            error "Please provide --fasta (the primary/target genome) when --perform_bbsplit is true."
        }

        // Static, user-curated list — read directly with Groovy rather
        // than through a Nextflow channel, since it's the same for
        // every sample and known up front, not per-sample data.
        def static_refs = []
        if (params.bbsplit_fasta_list) {
            file(params.bbsplit_fasta_list, checkIfExists: true).eachLine { line ->
                def (name, fasta_path) = line.tokenize(',')
                static_refs << [ name, file(fasta_path, checkIfExists: true) ]
            }
        }

        // Resolve the run-wide union of Kraken2-detected taxids to
        // reference fastas, merge with the static list, and shape it
        // into BBSplit's [other_ref_names, other_ref_paths] input.
        // toList() (not collect()!) guarantees exactly one emission —
        // an empty list [] if Kraken2 found nothing, or is off — so
        // this always pairs correctly against every sample in
        // ch_trimmed_reads. collect() looks equivalent but silently
        // emits NOTHING on an empty source channel, which would
        // starve BBMAP_BBSPLIT of its 4th input and it would never
        // run at all — with no error to say why.
        ch_bbsplit_refs = ch_kraken2_taxids
            .unique()
            .toList()
            .view { taxids -> "[bbsplit] unique taxid list collapsed via toList(): ${taxids}" } // TODO: Optional view for learning, marked for removal.
            .map { taxids ->
                def kraken2_refs = []
                if (taxids) {
                    if (!params.kraken2_taxid_lookup) {
                        error "Kraken2 detected contaminant taxa but --kraken2_taxid_lookup (CSV of taxid,name,fasta_path) was not provided; cannot resolve taxids to reference genomes."
                    }
                    def lookup = [:]
                    file(params.kraken2_taxid_lookup, checkIfExists: true).eachLine { line ->
                        def (taxid, name, fasta_path) = line.tokenize(',')
                        lookup[taxid] = [ name, fasta_path ]
                    }
                    taxids.each { taxid ->
                        if (lookup.containsKey(taxid)) {
                            def (name, fasta_path) = lookup[taxid]
                            kraken2_refs << [ name, file(fasta_path, checkIfExists: true) ]
                        } else {
                            log.warn "Kraken2 detected taxid ${taxid} but it has no entry in --kraken2_taxid_lookup — skipping it for BBSplit."
                        }
                    }
                }
                def all_refs = static_refs + kraken2_refs
                if (all_refs.isEmpty()) {
                    error "No contaminant references to screen against: --bbsplit_fasta_list is empty/unset and Kraken2 found no species-level taxa clearing the abundance/read-count thresholds. BBSplit needs at least one non-primary reference — either provide --bbsplit_fasta_list, or lower --kraken2_min_rel_abundance / --kraken2_min_reads, or disable --perform_bbsplit."
                }
                [ all_refs.collect { ref -> ref[0] }, all_refs.collect { ref -> ref[1] } ]
            }
            .view { refs -> "[bbsplit] final other_ref_names/other_ref_paths passed to BBMAP_BBSPLIT: ${refs}" } // TODO: Optional view for learning, marked for removal.

        BBMAP_BBSPLIT(
            ch_trimmed_reads,
            [], // index: none pre-built, build on-the-fly from the fastas below
            file(params.fasta, checkIfExists: true), // primary_ref: the human genome
            ch_bbsplit_refs, // other_ref_names, other_ref_paths
            false, // only_build_index: also split reads, not just build the index
        )

        ch_filtered_reads = BBMAP_BBSPLIT.out.primary_fastq
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBSPLIT.out.stats.map{ _meta, file -> file })
    }

    // ── STEP 5: SortMeRNA ──────────────────────────────────────
    // Ribosomal RNA depletion. 

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'rna-tools_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'rna-tools'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
