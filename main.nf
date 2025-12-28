#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = null          // S3 glob pattern or samplesheet CSV
params.outdir = null         // S3 output bucket/path
params.threads = 3
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

    // Stage VCF with its index file (.csi)
    vcf_with_index_ch = vcf_ch.map { vcf ->
         tuple(vcf, file(vcf.toUriString() + ".csi"))
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
    start_time=\$(date +%s.%N)
    Rscript ${rscript} ${gds_name} ${params.threads} ${variant_count} ${vcf}
    end_time=\$(date +%s.%N)
    awk -v start="\$start_time" -v end="\$end_time" -v variants="${variant_count}" \
        'BEGIN { printf "Conversion speed: %.2f variants/sec\\n", variants / (end - start) }'
    """
}
