#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {RUN_MINIMAP2; RUN_MMBWA_FROM_MINIMAP2; RUN_NGMLR; RUN_LAST; RUN_VACMAP} from './modules/aligners.nf'
include {RUN_SNIFFLES; RUN_CUTESV; RUN_DYSGU; RUN_DELLY; RUN_SAWFISH} from './modules/sv_callers.nf'
include {BENCHMARK} from './modules/benchmarking.nf'
include {PLOT_RESULTS} from './modules/plotting.nf'

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
params.regions_bed  = "/ref/chm13v2.0_censat_v2.1.bed"
params.dysgu_mode   = "nanopore-r10"
params.preset       = 'map-ont'
params.sample       = 'HG002'
params.ngmlr_mode   = 'ont'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    reads_ch = Channel.fromPath(params.reads)
    ref_ch = Channel.fromPath(params.ref)
    threads_ch = Channel.value(params.threads)
    regions_bed_ch = Channel.fromPath(params.regions_bed)
    sample_name = Channel.value(params.sample)
    preset_ch = Channel.value(params.preset)
    sq_lines_ch = Channel.fromPath(params.sq_lines)
    ngmlr_mode_ch = Channel.value(params.ngmlr_mode)

    RUN_MINIMAP2(ref_ch, reads_ch, threads_ch, sample_name, preset_ch)

    aln_ch = RUN_MINIMAP2.out.map { tuple -> tuple[2] }
    RUN_MMBWA_FROM_MINIMAP2(ref_ch, aln_ch, threads_ch, sample_name, regions_bed_ch, preset_ch)
    RUN_NGMLR(ref_ch, reads_ch, threads_ch, sample_name, ngmlr_mode_ch)
    RUN_LAST(ref_ch, reads_ch, threads_ch, sample_name, sq_lines_ch)
    RUN_VACMAP(ref_ch, reads_ch, threads_ch, sample_name)

    alignments = RUN_MINIMAP2.out
                .mix(RUN_MMBWA_FROM_MINIMAP2.out, RUN_NGMLR.out, RUN_LAST.out, RUN_VACMAP.out)

    RUN_SNIFFLES(ref_ch, alignments, params.platform)
    RUN_CUTESV(ref_ch, alignments, params.platform)
    RUN_DYSGU(ref_ch, alignments, params.platform, params.dysgu_mode)
    RUN_DELLY(ref_ch, alignments, params.platform)
    RUN_SAWFISH(ref_ch, alignments, params.platform)

    sv_callers_out = RUN_SNIFFLES.out.mix(
                     RUN_DYSGU.out,
                     RUN_CUTESV.out,
                     RUN_DELLY.out,
                     RUN_SAWFISH.out) // sample, aligner, platform, caller, vcf

    BENCHMARK(sv_callers_out, ref_ch, params.truth_vcf, params.truth_bed, params.bench_params)

    all_benchmarks = BENCHMARK.out.collect()

    PLOT_RESULTS(file("./plot_benchmark.py"), all_benchmarks)
}