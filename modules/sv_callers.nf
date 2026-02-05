/*
SV CALLERS
*/
process RUN_SNIFFLES {

    publishDir "${params.outdir}/calls/sniffles", mode: 'copy'
    tag "${sample}_${aligner}_sniffles"
    cpus { params.threads }

    input:
        path(ref)
        tuple val(sample), val(aligner), path(bam), path(bai)
        val(platform)

    output:
        tuple val(sample), val(aligner), val(platform), val('sniffles'), path("${sample}.${platform}.${aligner}.sniffles.vcf")

    script:
    """
    sniffles --threads ${task.cpus} --input ${bam} --vcf ${sample}.${platform}.${aligner}.sniffles.vcf
    """
}



process RUN_CUTESV {
    tag "${aligner}_cutesv"
    publishDir "${params.outdir}/calls/cutesv", mode: 'copy'
    cpus { params.threads }

    input:
        path(ref)
        tuple val(sample), val(aligner), path(bam), path(bai)
        val(platform)

    output:
        tuple val(sample), val(aligner), val(platform), val('cutesv'), path("${sample}.$platform.${aligner}.cutesv.vcf")

    script:
    """
    mkdir -p wd_cutesv
    cuteSV \
        --threads ${task.cpus} \
        -s 3 \
        ${bam} \
        ${ref} \
        ${sample}.$platform.${aligner}.cutesv.vcf \
        wd_cutesv
    """
}

process RUN_DYSGU {

    tag "${sample}_${aligner}_dysgu"
    publishDir "${params.outdir}/calls/dysgu", mode: 'copy'
    cpus { params.threads }

    input:
        path(ref)
        tuple val(sample), val(aligner), path(bam), path(bai)
        val(platform)
        val(dysgu_mode)

    output:
        tuple val(sample), val(aligner), val(platform), val('dysgu'), path("${sample}.${platform}.${aligner}.dysgu.vcf")

    script:
    """
    dysgu call \
        --mode ${dysgu_mode} \
        --procs ${task.cpus} \
        -x --clean \
        ${ref} wd ${bam} \
        > ${sample}.${platform}.${aligner}.dysgu.vcf
    """
}

process RUN_DELLY {

    tag "${sample}_${aligner}_delly"
    publishDir "${params.outdir}/calls/delly", mode: 'copy'
    cpus { params.threads }

    input:
    path(ref)
    tuple val(sample), val(aligner), path(bam), path(bai)
    val(platform)

    output:
        tuple val(sample), val(aligner), val(platform), val('delly'), path("${sample}.${platform}.${aligner}.delly.vcf")

    script:
    """
    delly lr -g ${ref} \
        ${bam} > ${sample}.${platform}.${aligner}.delly.vcf
    """
}

process RUN_SAWFISH {

    tag "${sample}_sawfish"
    publishDir "${params.outdir}/calls/sawfish", mode: 'copy'
    cpus { params.threads }

    input:
    path(ref)
    tuple val(sample), val(aligner), path(bam), path(bai)
    val(platform)

    when:
        platform == 'pacbio'

    output:
        tuple val(sample), val(aligner), val(platform), val('sawfish'), path("${sample}.${platform}.${aligner}.sawfish.vcf")

    script:
    """
    sawfish discover \
        --threads ${task.cpus} \
        --ref ${ref} \
        --bam ${bam} \
        --output-dir sawfish_discover

    sawfish joint-call \
        --threads ${task.cpus} \
        --sample sawfish_discover \
        --output-dir sawfish_call

    gunzip -c sawfish_call/genotyped.sv.vcf.gz \
        > ${sample}.${platform}.${aligner}.sawfish.vcf
    """
}
