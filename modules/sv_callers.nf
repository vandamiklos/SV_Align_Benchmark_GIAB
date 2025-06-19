process RUN_SNIFFLES {
    tag "${data_type}_${aligner}_sniffles"
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
    path ref
    path input_bam
    path idx
    val versions
    val data_type

    output:
    path "HG002.${data_type}.${aligner}.sniffles_*.vcf", emit: vcf

    script:
    def prefix = "HG002.${data_type}.${aligner}"
    """
    sniffles -t ${task.cpus} --input ${input_bam} --vcf ${prefix}.sniffles_${versions.sniffles}.vcf
    """
}
