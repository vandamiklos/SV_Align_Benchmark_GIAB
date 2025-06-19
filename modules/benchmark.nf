process BENCHMARK {
    publishDir "${params.outdir}/benchmarks", mode: 'copy'

    input:
    path ref
    path sv_call_vcfs
    path truth_vcf
    path truth_idx
    path truth_bed
    val versions
    val bench_params

    output:
    path "truv_output/*"

    script:
    """
    mkdir -p truv_output
    for vcf in ${sv_call_vcfs.join(' ')}; do
        truvari bench -f ${ref} -b ${truth_vcf} -c \$vcf -o truv_output/$(basename \$vcf .vcf) ${bench_params}
    done
    """
}
