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
        minimap2 -ax $preset -t $threads $ref $reads | samtools view -bh - | samtools sort -o '${sample}.minimap2.bam'
        samtools index '${sample}.minimap2.bam'
        """

}

process RUN_LAST {

    publishDir "${params.outdir}/alignment/last", mode: 'copy'

    tag "${sample}_last"
    cpus { threads }

    input:
        path ref
        path reads
        val threads
        val sample
        path sq_lines

    output:
        tuple val(sample), val('last'), path("${sample}.last.bam"), path("${sample}.last.bam.bai")

    script:
    """
    lastal --split -P $threads -Q0 $ref $reads | \
    maf-convert sam > '${sample}.lastalsplit.sam'; \
    samtools view -H '${sample}.lastalsplit.sam' > '${sample}.old_header.sam'; \
    cat '${sample}.old_header.sam' $sq_lines > '${sample}.new_header.sam'; \
    cat '${sample}.new_header.sam' '${sample}.lastalsplit.sam' | \
    samtools view -bh - | \
    samtools sort -o '${sample}.lastalsplit.bam'; \
    samtools index '${sample}.lastalsplit.bam'
    """
}

process RUN_NGMLR {

    publishDir "${params.outdir}/alignment/ngmlr", mode: 'copy'

    tag "${sample}_ngmlr"
    cpus { threads }

    input:
        path ref
        path reads
        val threads
        val sample
        val ngmlr_mode

    output:
        tuple val(sample), val('ngmlr'), path("${sample}.ngmlr.bam"), path("${sample}.ngmlr.bam.bai")

    script:
    """
    ngmlr -t $threads -r $ref -q $reads -x $ngmlr_mode | samtools view -bh - | samtools sort -o '${sample}.ngmlr.bam'
    samtools index '${sample}.ngmlr.bam'
    """
}

process RUN_VACMAP {

    publishDir "${params.outdir}/alignment/vacmap", mode: 'copy'

    tag "${sample}_vacmap"
    cpus { threads }

    input:
        path ref
        path reads
        val threads
        val sample

    output:
        tuple val(sample), val('vacmap'), path("${sample}.vacmap.bam"), path("${sample}.vacmap.bam.bai")

    script:
    """
    vacmap -ref $ref -read $reads -mode S -t $threads | samtools sort -o '${sample}.vacmap.bam'
    samtools index '${sample}.vacmap.bam'
    """
}
//
//process RUN_MMBWA {
//
//    publishDir "${params.outdir}/alignment/mmbwa", mode: 'copy'
//
//    tag "${sample}_mmbwa"
//   cpus { threads }
//
//    input:
//        path ref
//        path reads
//        val threads
//        val sample
//        val preset
//
//    output:
//        tuple val(sample), val('mmbwa'), path("${sample}.mmbwa.bam"), path("${sample}.mmbwa.bam.bai")
//
//    script:
//        """
//        mmbwa --input-fq $reads -t $threads --mm-args $preset --output mmbwa_out $ref
//
//        samtools view -bh mmbwa_out/final.bam |
//            samtools sort -o '${sample}.mmbwa.bam'
//
//        samtools index '${sample}.mmbwa.bam'
//       """
//}

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
        mmbwa --input-aln $aln -t $threads --mm-args $preset --regions_bed $bed_regions --output mmbwa_out $ref
        samtools view -bh mmbwa_out/final.bam | samtools sort -o '${sample}.mmbwa.bam'
        samtools index '${sample}.mmbwa.bam'
        """
}