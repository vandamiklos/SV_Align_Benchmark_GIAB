#!/usr/bin/env nextflow

process PLOT_RESULTS {
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path script
    path benchmark_dir

    output:
    path "plots/*"

    script:
    """
    python3 ${script} --input ${benchmark_dir} --output plots
    """
}
