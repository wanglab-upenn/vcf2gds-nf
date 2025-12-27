#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = null          // S3 glob pattern or samplesheet CSV
params.outdir = null         // S3 output bucket/path
params.threads = 6
params.container = "${HOME}/containers/r-base-gds-fuse-digest.4.5.1.sif"
params.script = "${moduleDir}/wl_script_R5.varcnt.R"

// Validate required parameters
if (!params.input) {
    error "Please specify input with --input (S3 glob pattern or samplesheet CSV)"
}
if (!params.outdir) {
    error "Please specify output directory with --outdir (e.g., s3://bucket/path)"
}

workflow {
    // Support both glob patterns and samplesheet CSV
    if (params.input.endsWith('.csv')) {
        vcf_ch = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.vcf, checkIfExists: true) }
    } else {
        vcf_ch = Channel.fromPath(params.input, checkIfExists: true)
    }

    // Stage VCF with its index file (.tbi or .csi)
    vcf_with_index_ch = vcf_ch.map { vcf ->
        def tbi = file("${vcf}.tbi")
        def csi = file("${vcf}.csi")
        def index = tbi.exists() ? tbi : (csi.exists() ? csi : null)
        return index ? tuple(vcf, index) : tuple(vcf, [])
    }

    script_ch = Channel.fromPath(params.script, checkIfExists: true)

    COUNT_VARIANTS(vcf_with_index_ch)
    CONVERT_TO_GDS(COUNT_VARIANTS.out, script_ch.collect())
}

process COUNT_VARIANTS {
    tag "${vcf.simpleName}"

    cpus 1
    memory '2G'
    // Load bcftools module before running
    beforeScript 'module load bcftools'

    input:
    tuple path(vcf), path(index)

    output:
    tuple path(vcf), path(index), env(variant_count)

    script:
    """
    # Index if needed (handles both .vcf.gz and .vcf.bgz)
    if [[ ! -f ${vcf}.csi && ! -f ${vcf}.tbi ]]; then
        bcftools index ${vcf}
    fi
    variant_count=\$(bcftools index --nrecords ${vcf})
    """
}

process CONVERT_TO_GDS {
    tag "${vcf.simpleName}"

    beforeScript 'module load apptainer'

    cpus params.threads
    memory '6G'

    container params.container

    publishDir "results", mode: 'copy'

    input:
    tuple path(vcf), path(index), val(variant_count)
    path rscript

    output:
    path "${gds_name}"

    script:
    // Handle both .vcf.gz and .vcf.bgz extensions
    gds_name = vcf.name.replaceAll(/\.vcf\.(gz|bgz)$/, '.vcf.gds')
    """
    Rscript ${rscript} ${gds_name} ${params.threads} ${variant_count} ${vcf}
    """
}
