#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INCLUDE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SETUP } from './modules/setup.nf'
include { RUN_BWA
          RUN_MINIMAP2
          RUN_LAST
          RUN_NGMLR
          RUN_VACMAP
          RUN_MMBWA
          RUN_MMBWA_FROM_MINIMAP2
 } from './modules/aligners.nf'
include {
    RUN_SNIFFLES
    RUN_CUTESV
    RUN_DYSGU
    RUN_DELLY
    RUN_SAWFISH
} from './modules/sv_callers.nf'
include { BENCHMARK }    from './modules/benchmark.nf'
include { PLOT_RESULTS } from './modules/plot.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.platform     = "ont"
params.threads      = 16
params.outdir       = "results"
params.input_reads  = "reads/*.fq"
params.ref          = "ref/chm13v2.0.fa"
params.truth_vcf    = "truth/HG002.SVs.vcf.gz"
params.truth_bed    = "truth/HG002.bed"
params.bench_params = "--passonly --giabreport --includebed ${params.truth_bed}"
params.sq_lines     = "ref/sq_lines.txt"
params.aligners     = ['bwa', 'minimap2', 'last', 'ngmlr', 'mmbwa']
params.sv_callers   = ['sniffles', 'cutesv', 'dysgu', 'delly', 'sawfish']

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

Channel.fromPath(params.input_reads).set { read_files }
Channel.value(file(params.ref)).set { ref }
Channel.value(params.threads).set { threads }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SETUP_WORKFLOW {
    main:
        SETUP()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGNMENT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

#workflow ALIGNMENT {
#
#    take:
#        sample
#        ref
#        reads
#        sq_lines
#        threads
#
#    main:
#        Channel.from(params.aligners).cross(sample).set {sample_aligner}
#
#        sample_aligner.filter { aligner, s -> aligner == 'bwa' }
#                      .map    { aligner, s -> s } | RUN_BWA(ref, reads, threads)
#
#        sample_aligner.filter { aligner, s -> aligner == 'minimap2' }
#                      .map    { aligner, s -> s } | RUN_MINIMAP2(ref, reads, threads)
#
#        sample_aligner.filter { aligner, s -> aligner == 'last' }
#                      .map    { aligner, s -> s } | RUN_LAST(ref, reads, sq_lines, threads)
#
#        sample_aligner.filter { aligner, s -> aligner == 'ngmlr' }
#                      .map    { aligner, s -> s } | RUN_NGMLR(ref, reads, threads)
#
#        sample_aligner.filter { aligner, s -> aligner == 'vacmap' }
#                      .map    { aligner, s -> s } | RUN_VACMAP(ref, reads, threads)
#
#        sample_aligner.filter { aligner, s -> aligner == 'mmbwa' }
#                      .map    { aligner, s -> s } | RUN_MMBWA(ref, reads, threads)
#
#    emit:
#        alignments =
#            RUN_BWA.out      \
#            .mix(RUN_MINIMAP2.out)
#            .mix(RUN_LAST.out)
#            .mix(RUN_NGMLR.out)
#            .mix(RUN_VACMAP.out)
#            .mix(RUN_MMBWA.out)
#}


workflow ALIGNMENT {

    take:
        sample
        ref
        reads
        threads

    main:
        minimap2_out = RUN_MINIMAP2(ref, reads, threads, sample)

        mmbwa_from_mm2 =
            minimap2_out
                .map { sample, aligner, bam, bai ->
                    tuple(ref, bam, threads, sample)
                }
                | RUN_MMBWA_FROM_MINIMAP2

        emit:
        alignments =
            minimap2_out
                .mix(mmbwa_from_mm2)
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SV CALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SV_CALLING {

    take:
        ref
        alignments   // (sample, aligner, bam, bai)
        versions = SETUP_WORKFLOW.out.versions

    main:
        sv_callers = params.sv_callers ?: []
        base =
            alignments
                .map { tuple(sample, aligner, bam, bai) ->
                    tuple(ref, bam, bai, versions, sample, aligner)
                }
        branches =
            base.branch {
                sniffles:sv_callers.contains('sniffles')
                cutesv:sv_callers.contains('cutesv')
                dysgu:sv_callers.contains('dysgu')
                delly:sv_callers.contains('delly') && it[5] != 'ngmlr'
                sawfish:sv_callers.contains('sawfish') && it[5] in ['minimap2','mmbwa']
            }

        if( sv_callers.contains('sniffles') )
            RUN_SNIFFLES(branches.sniffles)

        if( sv_callers.contains('cutesv') )
            RUN_CUTESV(branches.cutesv)

        if( sv_callers.contains('dysgu') )
            RUN_DYSGU(branches.dysgu)

        if( sv_callers.contains('delly') )
            RUN_DELLY(branches.delly)

        if( sv_callers.contains('sawfish') )
            RUN_SAWFISH(branches.sawfish)

    emit:
        sniffles = sv_callers.contains('sniffles') ? RUN_SNIFFLES.out.vcf : Channel.empty()
        cutesv   = sv_callers.contains('cutesv')   ? RUN_CUTESV.out.vcf   : Channel.empty()
        dysgu    = sv_callers.contains('dysgu')    ? RUN_DYSGU.out.vcf    : Channel.empty()
        delly    = sv_callers.contains('delly')    ? RUN_DELLY.out.vcf    : Channel.empty()
        sawfish  = sv_callers.contains('sawfish')  ? RUN_SAWFISH.out.vcf  : Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BENCHMARKING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BENCHMARKING {

    take:
        ref
        sv_calls     // (sample, aligner, caller, vcf)

    main:
        BENCHMARK(
            sv_calls,
            ref,
            file(params.truth_vcf),
            file(params.truth_vcf + ".tbi"),
            file(params.truth_bed),
            SETUP_WORKFLOW.out.versions,
            params.bench_params
        )
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PLOTTING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

    SETUP_WORKFLOW()

    ALIGNMENT(
        ref: ref,
        reads: read_files,
        threads: threads
    )

    SV_CALLING(
        ref: ref,
        alignments: ALIGNMENT.out
    )

    BENCHMARKING(
        ref: ref,
        sv_calls: SV_CALLING.out
    )

    PLOTTING(
        benchmark_dir: file("${params.outdir}/benchmarks")
    )
}

