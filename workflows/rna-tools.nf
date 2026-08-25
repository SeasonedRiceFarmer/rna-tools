/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
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
    // TODO: Might add another FASTQC step to ensure quality after FASTP. 

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
