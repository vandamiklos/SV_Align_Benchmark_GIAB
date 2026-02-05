/*
BENCHMARK
*/
process BENCHMARK {

    publishDir "${params.outdir}/benchmarks", mode: 'copy'

    input:
    tuple val(sample), val(aligner), val(platform), val(caller), path(vcf)
    path ref
    path truth_vcf
    path truth_bed
    val bench_params

    output:
    path("${sample}_truvari_GIAB_${platform}_${aligner}_${caller}")

    script:
    """
    bgzip -f ${vcf}
    tabix -f ${vcf}.gz
    mkdir -p ${sample}_truvari_GIAB_${platform}_${aligner}_${caller}
    truvari bench -f ${ref} -b ${truth_vcf} -c ${vcf}.gz -o ${sample}_truvari_GIAB_${platform}_${aligner}_${caller} ${bench_params} --includebed ${truth_bed}
    truvari refine -a wfa -f ${ref} --regions ${sample}_truvari_GIAB_${platform}_${aligner}_${caller}/candidate.refine.bed ${sample}_truvari_GIAB_${platform}_${aligner}_${caller}
    """
}