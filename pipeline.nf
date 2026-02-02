#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.platform     = "ont"
params.threads      = 16
params.outdir       = "/results"
params.reads        = "/reads/*.fq.gz"
params.ref          = "/ref/chm13v2.0.fa"
params.truth_vcf    = "/truth/GRCh38_HG2-T2TQ100-V1.1.vcf.gz"
params.truth_bed    = "/truth/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"
params.bench_params = "--passonly -r 1000 --dup-to-ins -p 0 --sizemax 50000"
params.sq_lines     = "/ref/sq_lines.txt"
params.aligners     = ['bwa', 'minimap2', 'last', 'ngmlr', 'mmbwa']
params.sv_callers   = ['sniffles', 'cutesv', 'dysgu', 'delly', 'sawfish']
params.regions_bed  = "/ref/chm13v2.0_censat_v2.1.bed"
params.dysgu_mode   = "nanopore-r10"
params.preset = 'map-ont'
params.sample = 'HG002'
/*
ALIGNERS
*/
//process RUN_BWA {
//
//   publishDir "${params.outdir}/alignment/bwa", mode: 'copy'
//
//    tag "${sample}_bwa"
//    cpus { threads }
//
//    input:
//        path ref
//        path reads
//        val threads
//        val sample
//
//    output:
//        tuple val(sample), val('bwa'), path("${sample}.bwa.bam"), path("${sample}.bwa.bam.bai")
//
//    script:
//    """
//    bwa mem -t $threads $ref $reads |
//        samtools view -bh - |
//        samtools sort -o '${sample}.bwa.bam'
//
//    samtools index '${sample}.bwa.bam'
//    """
//}

process RUN_MINIMAP2 {

    publishDir "${params.outdir}/alignment/minimap2", mode: 'copy'

    tag "${sample}_minimap2"
    cpus { threads }

    input:
        path ref
        path reads
        val threads
        val sample
        val preset

    output:
        tuple val(sample), val('minimap2'), path("${sample}.minimap2.bam"), path("${sample}.minimap2.bam.bai")
    script:
        """
        minimap2 -ax $preset -t $threads $ref $reads |
            samtools view -bh - |
            samtools sort -o '${sample}.minimap2.bam'

        samtools index '${sample}.minimap2.bam'
        """

}
//
//process RUN_LAST {
//
//    publishDir "${params.outdir}/alignment/last", mode: 'copy'
//
//    tag "${sample}_last"
//    cpus { threads }
//
//    input:
//        path ref
//        path reads
//        path sq_lines
//        val threads
//        val sample
//
//    output:
//        tuple val(sample), val('last'), path("${sample}.last.bam"), path("${sample}.last.bam.bai")
//
//    script:
//    """
//    head -n 100 $reads > training_set.fq
//    last-train -Q0 $ref training_set.fq > trained_parameters.txt
//
//    lastal --split -p trained_parameters.txt -P $threads -Q0 $ref $reads |
//        maf-convert sam > last.sam
//
//    samtools view -H last.sam > old_header.sam
//    cat old_header.sam $sq_lines > new_header.sam
//    cat new_header.sam last.sam |
//        samtools view -bh - |
//        samtools sort -o '${sample}.last.bam'
//
//    samtools index '${sample}.last.bam'
//    """
//}
//
//process RUN_NGMLR {
//
//    publishDir "${params.outdir}/alignment/ngmlr", mode: 'copy'
//
//    tag "${sample}_ngmlr"
//    cpus { threads }
//
//    input:
//        path ref
//        path reads
//        val threads
//        val sample
//
//    output:
//        tuple val(sample), val('ngmlr'), path("${sample}.ngmlr.bam"), path("${sample}.ngmlr.bam.bai")
//
//    script:
//    """
//    ngmlr -t $threads -r $ref -q $reads |
//        samtools view -bh - |
//        samtools sort -o '${sample}.ngmlr.bam'
//
//    samtools index '${sample}.ngmlr.bam'
//    """
//}
//
//process RUN_VACMAP {
//
//    publishDir "${params.outdir}/alignment/vacmap", mode: 'copy'
//
//    tag "${sample}_vacmap"
//    cpus { threads }
//
//    input:
//        path ref
//        path reads
//        val threads
//        val sample
//
//    output:
//        tuple val(sample), val('vacmap'), path("${sample}.vacmap.bam"), path("${sample}.vacmap.bam.bai")
//
//    script:
//    """
//    vacmap -ref $ref -read $reads -mode S -t $threads |
//        samtools sort -o '${sample}.vacmap.bam'
//
//    samtools index '${sample}.vacmap.bam'
//    """
//}
//
process RUN_MMBWA {

    publishDir "${params.outdir}/alignment/mmbwa", mode: 'copy'

    tag "${sample}_mmbwa"
    cpus { threads }

    input:
        path ref
        path reads
        val threads
        val sample
        val preset

    output:
        tuple val(sample), val('mmbwa'), path("${sample}.mmbwa.bam"), path("${sample}.mmbwa.bam.bai")

    script:
        """
        mmbwa --input-fq $reads -t $threads --mm-args $preset --output mmbwa_out $ref

        samtools view -bh mmbwa_out/final.bam |
            samtools sort -o '${sample}.mmbwa.bam'

        samtools index '${sample}.mmbwa.bam'
        """
}

process RUN_MMBWA_FROM_MINIMAP2 {

    publishDir "${params.outdir}/alignment/mmbwa", mode: 'copy'

    tag "${sample}_mmbwa"
    cpus { threads }

    input:
        path ref
        path aln
        val threads
        val sample
        path bed_regions
        val preset

    output:
        tuple val(sample), val('mmbwa'), path("${sample}.mmbwa.bam"), path("${sample}.mmbwa.bam.bai")

    script:
        """
        mmbwa --input-aln $aln \
              -t $threads \
              --mm-args $preset \
              --regions_bed $bed_regions \
              --output mmbwa_out $ref

        samtools view -bh mmbwa_out/final.bam |
            samtools sort -o '${sample}.mmbwa.bam'

        samtools index '${sample}.mmbwa.bam'
        """
}

/*
SV CALLERS
*/
process RUN_SNIFFLES {

    publishDir "${params.outdir}/calls/sniffles", mode: 'copy'
    tag "${sample}_${aligner}_sniffles"
    cpus { params.threads }

    input:
        tuple path(ref), path(bam), path(bai), val(sample), val(aligner), val(platform), val(dysgu_mode)

    output:
        tuple val(sample), val(aligner), val('sniffles'), path("${sample}.${platform}.${aligner}.sniffles.vcf")

    script:
    """
    sniffles --threads ${task.cpus} \
        --input ${bam} \
        --vcf ${sample}.${platform}.${aligner}.sniffles.vcf
    """
}



process RUN_CUTESV {
    tag "${aligner}_cutesv"
    publishDir "${params.outdir}/calls/cutesv", mode: 'copy'
    cpus { params.threads }

    input:
    tuple path(ref), path(bam), path(bai), val(sample), val(aligner), val(platform), val(dysgu_mode)

    output:
    tuple val(sample), val(aligner), val('cutesv'), path("${sample}.$platform.${aligner}.cutesv.vcf")

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
    tuple path(ref), path(bam), path(bai), val(sample), val(aligner), val(platform), val(dysgu_mode)

    output:
    tuple val(sample), val(aligner), val('dysgu'), path("${sample}.${platform}.${aligner}.dysgu.vcf")

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
    tuple path(ref), path(bam), path(bai), val(sample), val(aligner), val(platform), val(dysgu_mode)

    output:
    tuple val(sample), val(aligner), val('delly'), path("${sample}.${platform}.${aligner}.delly.vcf")

    script:
    """
    delly lr -g ${ref} \
        ${bam} \
        > ${sample}.${platform}.${aligner}.delly.vcf
    """
}

process RUN_SAWFISH {

    tag "${sample}_sawfish"
    publishDir "${params.outdir}/calls/sawfish", mode: 'copy'
    cpus { params.threads }

    input:
    tuple path(ref), path(bam), path(bai), val(sample), val(aligner), val(platform), val(dysgu_mode)

    when:
        platform == 'pacbio'

    output:
    tuple val(sample), val(aligner), val('sawfish'), path("${sample}.${platform}.${aligner}.sawfish.vcf")

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

/*
BENCHMARK
*/

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

/*
PLOT
*/
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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow ALIGNMENT {

    take:
        ref
        reads
        threads
        sample
        bed_regions

    main:
        minimap2_out =
            RUN_MINIMAP2(ref, reads, threads, sample)

        mmbwa_from_mm2 =
            minimap2_out
                .map { sample, aligner, bam, bai ->
                    tuple(ref, bam, threads, sample, bed_regions)
                }
                | RUN_MMBWA_FROM_MINIMAP2
                .map { sample, aligner, bam, bai ->
                    tuple(sample, aligner, bam, bai)
                }

    emit:
        alignments =
            minimap2_out
                .mix(mmbwa_from_mm2)
}


workflow SV_CALLING {

    take:
        ref
        alignments   // (sample, aligner, bam, bai)

    main:
        sv_callers = params.sv_callers ?: []

        base =
            alignments.map { sample, aligner, bam, bai ->
                [
                    ref: ref,
                    bam: bam,
                    bai: bai,
                    sample: sample,
                    aligner: aligner,
                    platform: params.platform,
                    dysgu_mode: params.dysgu_mode
                ]
            }

        /*
         * Aligner-specific filters
         */
        non_ngmlr   = base.filter { it.aligner != 'ngmlr' }
        mm2_mmbwa   = base.filter { it.aligner in ['minimap2', 'mmbwa'] }

        /*
         * SV callers
         */
        if( sv_callers.contains('sniffles') )
            RUN_SNIFFLES(base)

        if( sv_callers.contains('cutesv') )
            RUN_CUTESV(base)

        if( sv_callers.contains('dysgu') )
            RUN_DYSGU(base)

        if( sv_callers.contains('delly') )
            RUN_DELLY(non_ngmlr)

        if( sv_callers.contains('sawfish') )
            RUN_SAWFISH(mm2_mmbwa)

    emit:
        sniffles = sv_callers.contains('sniffles') ? RUN_SNIFFLES.out : Channel.empty()
        cutesv   = sv_callers.contains('cutesv')   ? RUN_CUTESV.out   : Channel.empty()
        dysgu    = sv_callers.contains('dysgu')    ? RUN_DYSGU.out    : Channel.empty()
        delly    = sv_callers.contains('delly')    ? RUN_DELLY.out    : Channel.empty()
        sawfish  = sv_callers.contains('sawfish')  ? RUN_SAWFISH.out  : Channel.empty()
}


workflow BENCHMARKING {

    take:
        ref
        sv_calls     // (sample, aligner, caller, vcf)
        platform

    main:
        BENCHMARK(
            sv_calls,
            ref,
            file(params.truth_vcf),
            file(params.truth_vcf + ".tbi"),
            file(params.truth_bed),
            params.bench_params,
            params.platform
        )
    emit:
        out = BENCHMARK.out
}

workflow PLOTTING {

    take:
        benchmark_dir

    main:
        PLOT_RESULTS(
            file("./plot_benchmark.py"),
            benchmark_dir
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FINAL PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    reads_ch = Channel.fromPath(params.reads)
    ref_ch = Channel.fromPath(params.ref)
    threads_ch = Channel.value(params.threads)
    regions_bed_ch = Channel.fromPath(params.regions_bed)
    sample_name = Channel.value(params.sample)
    preset_ch = Channel.value(params.preset)

    /*SETUP_WORKFLOW()*/
    RUN_MINIMAP2(ref_ch, reads_ch, threads_ch, sample_name, preset_ch)

    aln_ch = RUN_MINIMAP2.out.map { tuple -> tuple[2] }
    RUN_MMBWA_FROM_MINIMAP2(ref_ch, aln_ch, threads_ch, sample_name, regions_bed_ch, preset_ch)

    alignments = RUN_MINIMAP2.out
                .mix(RUN_MMBWA_FROM_MINIMAP2.out)

    SV_CALLING(
        ref: ref_ch,
        alignments: alignments
    )

    BENCHMARKING(
        ref: ref_ch,
        sv_calls: SV_CALLING.out
    )

    all_benchmarks = BENCHMARKING.out.collect()

    PLOTTING(
        benchmark_dir: all_benchmarks
    )
}