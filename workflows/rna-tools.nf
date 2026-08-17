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
    //
    // MODULE: Run FastQC
    //
    // ch_samplesheet is built upstream (in PIPELINE_INITIALISATION, from your
    // samplesheet.csv) and has the shape:
    //     [ meta, reads ]
    // where:
    //   meta  = a Groovy map, e.g. [ id:'WHOLE_BLOOD_TEST', single_end:false ]
    //   reads = a *list* of fastq.gz paths — one item if single-end, two if paired-end
    // Every process below re-uses this same [ meta, reads ] shape as its input,
    // and (almost) every process re-emits [ meta, <its own output files> ] so the
    // metadata (sample id, single_end, later: condition/replicate info) travels
    // alongside the data all the way through the pipeline instead of getting lost.
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // MODULE: Run fastp (adapter/quality trimming)
    //
    // FastQC only *looked* at the raw reads and wrote a report — it did not change
    // any files. fastp is the step that actually edits the reads: it cuts adapter
    // sequence off the ends, drops low-quality bases/reads, and writes new,
    // cleaned fastq.gz files. STAR (added in the next step) will align these
    // cleaned reads, not the raw ones.
    //
    // The FASTP module (modules/nf-core/fastp/main.nf) expects THREE things per
    // sample instead of two:
    //     tuple val(meta), path(reads), path(adapter_fasta)
    // adapter_fasta lets you supply a known list of adapter sequences to trim.
    // We pass an empty list [] for it, which tells fastp "auto-detect the adapter
    // instead" — fastp is good at this, and it's the standard default. So we
    // .map() over ch_samplesheet just to bolt an empty list onto each tuple.
    ch_fastp_input = ch_samplesheet.map { meta, reads -> [ meta, reads, [] ] }

    // The three plain `false` arguments below are NOT per-sample data — they're
    // single shared settings for every sample, matching how the module itself is
    // written (input[1], input[2], input[3] in its test file):
    //   discard_trimmed_pass : false = actually keep/emit the trimmed reads
    //                          (if this were true, fastp would trim but not save them)
    //   save_trimmed_fail    : false = don't bother keeping reads that failed QC
    //   save_merged          : false = don't merge overlapping paired-end mates
    //                          into one read (a feature for very short fragments;
    //                          not relevant to standard RNA-seq)
    FASTP(
        ch_fastp_input,
        false, // discard_trimmed_pass
        false, // save_trimmed_fail
        false, // save_merged
    )

    // FASTP.out.reads = [ meta, trimmed_reads ] — this is what STAR will consume next.
    // We don't use it yet this step, but we name it here so it's obvious where the
    // pipeline picks back up.
    ch_trimmed_reads = FASTP.out.reads

    // fastp writes its own QC report as a .json file. MultiQC knows how to read
    // fastp's json format natively, so we just add it to the same "pile of files"
    // channel (ch_multiqc_files) that FastQC's output already goes into. MultiQC
    // will scan every file in this pile at the end and figure out which tool made
    // each one — we don't have to tell it explicitly.
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
