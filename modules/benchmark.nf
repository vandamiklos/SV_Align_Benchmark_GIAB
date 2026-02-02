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
    val platform

    output:
    path("truvari_GIAB_${platform}_${caller}_${sample}")

    script:
    """
    bgzip -f \${vcf}
    tabix -f \${vcf}.gz
    out_dir=truvari_GIAB_${platform}_${caller}_${sample}
    mkdir -p $out_dir
    truvari bench \
        -f ${ref} \
        -b ${truth_vcf} \
        -c ${vcf}.gz \
        -o $out_dir \
        ${bench_params} \
        --includebed ${truth_bed}

    truvari refine -a wfa -f ${ref} \
                   --regions ${out_dir}/candidate.refine.bed \
                   ${out_dir}
    """
}
