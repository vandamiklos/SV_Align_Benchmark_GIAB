#!/usr/bin/env nextflow
process BENCHMARK {

    publishDir "${params.outdir}/benchmarks", mode: 'copy'

    input:
    tuple val(sample), val(aligner), val(caller), path(vcf)
    path ref
    path truth_vcf
    path truth_idx
    path truth_bed
    val versions
    val bench_params

    output:
    tuple val(sample), val(aligner), val(caller), path("truv_output")

    script:
    """
    mkdir -p truv_output
    truvari bench \
        -f ${ref} \
        -b ${truth_vcf} \
        -c ${vcf} \
        -o truv_output \
        ${bench_params}
    """
}
