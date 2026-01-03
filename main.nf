#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { validateParameters ; paramsSummaryLog ; samplesheetToList } from 'plugin/nf-schema'

// S3 output bucket/path
params.outdir = null

params.threads = 3
params.container = "${HOME}/images/r-base-gds-fuse-digest.4.5.1.sif"
params.script = "${moduleDir}/wl_script_R5.varcnt.R"

def validateParams() {
    // Check that either sample sheet or single sample is provided
    if (!params.input_csv && !params.sample) {
        error("ERROR: You must provide either --input_csv (sample sheet CSV) or --sample")
    }

    // Validate sample sheet exists if provided
    if (params.input_csv) {
        def samplesheet = file(params.input_csv)
        if (!samplesheet.exists()) {
            error("ERROR: Sample sheet not found: ${params.input_csv}")
        }
    }
}


workflow {

    log.info(paramsSummaryLog(workflow))

    // Validate parameters
    validateParams()

    // Validate required parameters
    // if (!params.input_csv) {
    //     error("Please specify input with --input_csv (S3 glob pattern or samplesheet CSV)")
    // }
    if (!params.outdir) {
        error("Please specify output directory with --outdir (e.g., s3://bucket/path)")
    }

    // Support both glob patterns and samplesheet CSV
    if (params.input_csv.endsWith('.csv')) {
        vcf_ch = channel.fromPath(params.input_csv)
            .splitCsv(header: false)
            .map { row -> file(row[0], checkIfExists: true) }
    }
    else {
        vcf_ch = channel.fromPath(params.input_s3, checkIfExists: true)
    }

    // Stage VCF with its index file (.csi)
    vcf_with_index_ch = vcf_ch.map { vcf ->
        tuple(vcf, file(vcf.toUriString() + ".csi"))
    }

    script_ch = channel.fromPath(params.script, checkIfExists: true)

    COUNT_VARIANTS(vcf_with_index_ch)
    CONVERT_TO_GDS(COUNT_VARIANTS.out, script_ch.collect())

    workflow.onComplete = {

        if (!workflow.success) {
            return null
        }

        // Clean up S3 staged files
        def stageDir = "${workflow.workDir}/stage-${workflow.sessionId}"
        log.info("Cleaning up staged files at: ${stageDir}")

//        try {
            def cmd = "rm -f --recursive ${stageDir}"
            def proc = cmd.execute()
            proc.waitFor()

            if (proc.exitValue() == 0) {
                log.info("Successfully cleaned up staged files")
            }  else {
                log.warn("Stage cleanup returned non-zero exit: ${proc.err.text}")
            }
//        } catch (e: Exception) {
//            log.warn("Failed to clean up staged files: ${e.message}")
//        }
    }// end oncomplete
}

process COUNT_VARIANTS {
    tag "${vcf.simpleName}"

    cpus 1
    memory '2G'

    executor 'local'

    input:
    tuple path(vcf), path(index)

    output:
    tuple path(vcf), path(index), env(variant_count)

    script:
    def vcf_name = vcf.name
    def gz_name = vcf_name.replaceAll(/\.bgz$/, '.gz')
    def is_bgz = vcf_name.endsWith('.bgz')
    """
    ls -alLtFh

    if [[ "${is_bgz}" == "true" ]]; then
        # Symlink .bgz to .gz (bcftools doesn't recognize .bgz)
        ln -sf ${vcf} "${gz_name}"
        ln -sf ${index} "${gz_name}.csi"
        variant_count=\$(bcftools index --nrecords "${gz_name}")
    else
        variant_count=\$(bcftools index --nrecords ${vcf})
    fi

    """
}

process CONVERT_TO_GDS {
    tag "${vcf.simpleName}"

    cpus params.threads
    memory '6G'

    container params.container

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(vcf), path(index), val(variant_count)
    path rscript

    output:
    path "${gds_name}"

    script:
    // Handle both .vcf.gz and .vcf.bgz extensions
    gds_name = vcf.name.replaceAll(/\.vcf\.(gz|bgz)$/, '.vcf.gds')
    """
    start_time=\$(date +%s.%N)
    Rscript ${rscript} ${gds_name} ${params.threads} ${variant_count} ${vcf}
    end_time=\$(date +%s.%N)
    awk -v start="\$start_time" -v end="\$end_time" -v variants="${variant_count}" \
        'BEGIN { printf "Conversion speed: %.2f variants/sec\\n", variants / (end - start) }'
    """
}
