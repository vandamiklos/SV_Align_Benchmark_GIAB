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
        val aligner
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'sniffles',
        path "${sample}.${aligner}.sniffles.vcf"
    ),
    emit: vcf

    script:
    """
    sniffles \
        --threads ${task.cpus} \
        --input ${bam} \
        --vcf ${sample}.${aligner}.sniffles.vcf
    """
}

process RUN_CUTESV {

    tag "${sample}_${aligner}_cutesv"
    publishDir "${params.outdir}/calls/cutesv", mode: 'copy'
    cpus { params.threads }

    input:
    tuple(
        path ref,
        path bam,
        path bai,
        val sample,
        val aligner
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'cutesv',
        path "${sample}.${aligner}.cutesv.vcf"
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
        ${sample}.${aligner}.cutesv.vcf \
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
        val aligner
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'dysgu',
        path "${sample}.${aligner}.dysgu.vcf"
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
        > ${sample}.${aligner}.dysgu.vcf
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
        val aligner
    )

    output:
    tuple(
        val sample,
        val aligner,
        val 'delly',
        path "${sample}.${aligner}.delly.vcf"
    ),
    emit: vcf

    script:
    def tech = params.caller_params[aligner]?.delly?.tech ?: 'ONT'
    """
    delly lr \
        --technology ${tech} \
        -g ${ref} \
        ${bam} \
        > ${sample}.${aligner}.delly.vcf
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
        val aligner
    )

    when:
    aligner == 'pacbio'

    output:
    tuple(
        val sample,
        val aligner,
        val 'sawfish',
        path "${sample}.${aligner}.sawfish.vcf"
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
        > ${sample}.${aligner}.sawfish.vcf
    """
}

