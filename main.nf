nextflow.enable.dsl=2

// Include modules
include { SETUP } from './modules/setup.nf'
include { AlignReads } from './modules/aligners.nf'
include {
    RUN_SNIFFLES;
    RUN_CUTESV;
    RUN_DYSGU;
    RUN_DELLY;
    RUN_SAWFISH
} from './modules/sv_callers.nf'
include { BENCHMARK } from './modules/benchmark.nf'
include { PLOT_RESULTS } from './modules/plot.nf'

// Input params
params.platform = "ont"
params.threads = 22
params.outdir = "results"
params.input_reads = "reads/*.fq"
params.ref = "ref/chm13v2.0.fa"
params.truth_vcf = "truth/HG002.SVs.vcf.gz"
params.truth_bed = "truth/HG002.bed"
params.bench_params = "--passonly --giabreport --includebed ${params.truth_bed}"

// Channels
Channel.fromPath(params.input_reads).set { read_files }
Channel.value(params.ref).set { ref }
Channel.value(params.threads).set { threads }
Channel.value(params.outdir).set { outdir }

// Setup step
workflow SETUP_WORKFLOW {
    main:
        SETUP()
}

// Align reads with all aligners
workflow ALIGNMENT {
    take:
        ref
        reads = read_files
        threads
        outdir

    main:
        Channel.of("bwa", "minimap2", "last", "ngmlr", "vacmap").set { aligners }
        aligners
            .combine(ref)
            .combine(reads)
            .combine(SETUP_WORKFLOW.out.sq_lines)
            .combine(threads)
            .combine(outdir)
            | map { it.flatten() }
            | AlignReads
}

// SV calling
workflow SV_CALLING {
    take:
        ref
        alignments
        versions = SETUP_WORKFLOW.out.versions

    main:
        alignments
            .map { bam_file ->
                def aligner = bam_file.baseName.split(/\./)[1] // e.g., HG002.ont.bwa.bam
                def data_type = bam_file.baseName.split(/\./)[0] // e.g., HG002.ont.bwa.bam
                tuple(ref, bam_file, "${bam_file}.bai", versions, data_type, aligner)
            }
            | flatten()
            | branch {
                sniffles: it | RUN_SNIFFLES
                cutesv:   it | RUN_CUTESV
                dysgu:    it | RUN_DYSGU
                delly:    it | RUN_DELLY
                sawfish:  it | RUN_SAWFISH
            }
}

// Benchmarking
workflow BENCHMARKING {
    take:
        ref
        sv_calls
        truth_vcf = file(params.truth_vcf)
        truth_idx = file(params.truth_vcf + ".tbi")
        truth_bed = file(params.truth_bed)

    main:
        BENCHMARK(ref, sv_calls.flatten(), truth_vcf, truth_idx, truth_bed, SETUP_WORKFLOW.out.versions, params.bench_params)
}

// Plotting
workflow PLOTTING {
    take:
        benchmark_dir = BENCHMARKING.out

    main:
        PLOT_RESULTS(file("scripts/plot_truvari.py"), benchmark_dir)
}

// Final pipeline
workflow {
    SETUP_WORKFLOW()
    ALIGNMENT(ref: ref, reads: read_files, threads: threads, outdir: outdir)
    SV_CALLING(ref: ref, alignments: ALIGNMENT.out)
    BENCHMARKING(ref: ref, sv_calls: SV_CALLING.out)
    PLOTTING(benchmark_dir: BENCHMARKING.out)
}

