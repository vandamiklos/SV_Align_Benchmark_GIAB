process BENCHMARK {

    publishDir "${params.outdir}/benchmarks", mode: 'copy'

    input:
    tuple val(sample), val(aligner), val(caller), path(vcf)
    path ref
    path truth_vcf
    path truth_idx
    path truth_bed
    val bench_params
    val platform

    output:
    path("truvari_GIAB_${platform}_${caller}_${sample}")

    script:
    """
    out_dir=truvari_GIAB_${platform}_${caller}_${sample}
    mkdir -p $out_dir
    truvari bench \
        -f ${ref} \
        -b ${truth_vcf} \
        -c ${vcf} \
        -o $out_dir \
        ${bench_params}
    """
}
