#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process RUN_SNIFFLES {

    tag "${sample}_${aligner}_sniffles"
    publishDir "${params.outdir}/calls/sniffles", mode: 'copy'

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val versions,
        val sample,
        val aligner
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'sniffles',
        path "${sample}.${aligner}.sniffles_${versions.sniffles}.vcf"
    ),
    emit: vcf

    script:
    """
    sniffles \
        --threads ${task.cpus} \
        --input ${bam} \
        --vcf ${sample}.${aligner}.sniffles_${versions.sniffles}.vcf
    """
}

process RUN_CUTESV {

    tag "${sample}_${aligner}_cutesv"
    publishDir "${params.outdir}/calls/cutesv", mode: 'copy'

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val versions,
        val sample,
        val aligner
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'cutesv',
        path "${sample}.${aligner}.cutesv_${versions.cutesv}.vcf"
    ),
    emit: vcf

    script:
    """
    mkdir -p wd_cutesv
    cuteSV \
        --threads ${task.cpus} \
        ${params.cutesv_opts ?: ''} \
        ${bam} \
        ${ref} \
        ${sample}.${aligner}.cutesv_${versions.cutesv}.vcf \
        wd_cutesv
    """
}

process RUN_DYSGU {

    tag "${sample}_${aligner}_dysgu"
    publishDir "${params.outdir}/calls/dysgu", mode: 'copy'

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val versions,
        val sample,
        val aligner
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'dysgu',
        path "${sample}.${aligner}.dysgu_${versions.dysgu}.vcf"
    ),
    emit: vcf

    script:
    def mode = params.caller_params[aligner]?.dysgu?.mode ?: 'pacbio'
    """
    dysgu call \
        --mode ${mode} \
        --procs ${task.cpus} \
        -x --clean \
        ${ref} wd ${bam} \
        > ${sample}.${aligner}.dysgu_${versions.dysgu}.vcf
    """
}

process RUN_DELLY {

    tag "${sample}_${aligner}_delly"
    publishDir "${params.outdir}/calls/delly", mode: 'copy'

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val versions,
        val sample,
        val aligner
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'delly',
        path "${sample}.${aligner}.delly_${versions.delly}.vcf"
    ),
    emit: vcf

    script:
    def tech = params.caller_params[aligner]?.delly?.tech ?: 'ONT'
    """
    delly lr \
        --technology ${tech} \
        -g ${ref} \
        ${bam} \
        > ${sample}.${aligner}.delly_${versions.delly}.vcf
    """
}

process RUN_SAWFISH {

    tag "${sample}_sawfish"
    publishDir "${params.outdir}/calls/sawfish", mode: 'copy'

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val versions,
        val sample,
        val aligner
    )

    when:
    aligner == 'pacbio'

    output:
    tuple(
        val sample,
        val aligner,
        val 'sawfish',
        path "${sample}.${aligner}.sawfish_${versions.sawfish}.vcf"
    ),
    emit: vcf

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
        > ${sample}.${aligner}.sawfish_${versions.sawfish}.vcf
    """
}

