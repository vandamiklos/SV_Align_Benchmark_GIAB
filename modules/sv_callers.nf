process RUN_SNIFFLES {

    tag "${sample}_${aligner}_sniffles"
    publishDir "${params.outdir}/calls/sniffles", mode: 'copy'
    cpus { params.threads }

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val sample,
        val aligner,
        val platform,
        val dysgu_mode
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'sniffles',
        path "${sample}.${platform}.${aligner}.sniffles.vcf"
    ),
    emit: vcf

    script:
    """
    sniffles \
        --threads ${task.cpus} \
        --input ${bam} \
        --vcf ${sample}.${platform}.${aligner}.sniffles.vcf
    """
}

process RUN_CUTESV {
    tag "${aligner}_cutesv"
    publishDir "${params.outdir}/calls/cutesv", mode: 'copy'
    cpus { params.threads }

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val sample,
        val aligner,
        val platform,
        val dysgu_mode
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'cutesv',
        path "${sample}.$platform.${aligner}.cutesv.vcf"
    ),
    emit: vcf

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
    tuple(
        path ref,
        path bam,
        path bai,
        val sample,
        val aligner,
        val platform,
        val dysgu_mode
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'dysgu',
        path "${sample}.${platform}.${aligner}.dysgu.vcf"
    ),
    emit: vcf

    script:
    def mode = params.platform.dysgu.mode
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
    tuple(
        path ref,
        path bam,
        path bai,
        val sample,
        val aligner,
        val platform,
        val dysgu_mode
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'delly',
        path "${sample}.${platform}.${aligner}.delly.vcf"
    ),
    emit: vcf

    script:
    def tech = params.caller_params[aligner]?.delly?.tech ?: 'ONT'
    """
    delly lr \
        --technology ${tech} \
        -g ${ref} \
        ${bam} \
        > ${sample}.${platform}.${aligner}.delly.vcf
    """
}

process RUN_SAWFISH {

    tag "${sample}_sawfish"
    publishDir "${params.outdir}/calls/sawfish", mode: 'copy'
    cpus { params.threads }

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val sample,
        val aligner,
        val platform,
        val dysgu_mode
    )

    when:
    aligner == 'pacbio'

    output:
    tuple(
        val sample,
        val aligner,
        val 'sawfish',
        path "${sample}.${platform}.${aligner}.sawfish.vcf"
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
        > ${sample}.${platform}.${aligner}.sawfish.vcf
    """
}

