process RUN_BWA {

    publishDir "${params.outdir}/alignment/bwa", mode: 'copy'

    tag "${sample}_bwa"
    cpus { threads }

    input:
        path ref
        path reads
        val threads
        val sample

    output:
        tuple val(sample), val('bwa'), path("${sample}.bwa.bam"), path("${sample}.bwa.bam.bai")

    script:
    """
    bwa mem -t $threads $ref $reads |
        samtools view -bh - |
        samtools sort -o '${sample}.bwa.bam'

    samtools index '${sample}.bwa.bam'
    """
}

process RUN_MINIMAP2 {

    publishDir "${params.outdir}/alignment/minimap2", mode: 'copy'

    tag "${sample}_minimap2"
    cpus { threads }

    input:
        path ref
        path reads
        val threads
        val sample

    output:
        tuple val(sample), val('minimap2'), path("${sample}.minimap2.bam"), path("${sample}.minimap2.bam.bai")

    script:
    def preset = params.platform == 'ont' ? 'map-ont' : 'map-hifi'
    """
    minimap2 -ax $preset -t $threads $ref $reads |
        samtools view -bh - |
        samtools sort -o '${sample}.minimap2.bam'

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
        path sq_lines
        val threads
        val sample

    output:
        tuple val(sample), val('last'), path("${sample}.last.bam"), path("${sample}.last.bam.bai")

    script:
    """
    head -n 100 $reads > training_set.fq
    last-train -Q0 $ref training_set.fq > trained_parameters.txt

    lastal --split -p trained_parameters.txt -P $threads -Q0 $ref $reads |
        maf-convert sam > last.sam

    samtools view -H last.sam > old_header.sam
    cat old_header.sam $sq_lines > new_header.sam
    cat new_header.sam last.sam |
        samtools view -bh - |
        samtools sort -o '${sample}.last.bam'

    samtools index '${sample}.last.bam'
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

    output:
        tuple val(sample), val('ngmlr'), path("${sample}.ngmlr.bam"), path("${sample}.ngmlr.bam.bai")

    script:
    """
    ngmlr -t $threads -r $ref -q $reads |
        samtools view -bh - |
        samtools sort -o '${sample}.ngmlr.bam'

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
    vacmap -ref $ref -read $reads -mode S -t $threads |
        samtools sort -o '${sample}.vacmap.bam'

    samtools index '${sample}.vacmap.bam'
    """
}

process RUN_MMBWA {

    publishDir "${params.outdir}/alignment/mmbwa", mode: 'copy'

    tag "${sample}_mmbwa"
    cpus { threads }

    input:
        path ref
        path reads
        val threads
        val sample

    output:
        tuple val(sample), val('mmbwa'), path("${sample}.mmbwa.bam"), path("${sample}.mmbwa.bam.bai")

    script:
    def preset = params.platform == 'ont' ? 'map-ont' : 'map-hifi'
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

    output:
        tuple val(sample), val('mmbwa'), path("${sample}.mmbwa.bam"), path("${sample}.mmbwa.bam.bai")

    script:
    def preset = params.platform == 'ont' ? 'map-ont' : 'map-hifi'
    """
    mmbwa --input-aln $aln -t $threads --mm-args $preset --bed_regions $bed_regions --output mmbwa_out $ref

    samtools view -bh mmbwa_out/final.bam |
        samtools sort -o '${sample}.mmbwa.bam'

    samtools index '${sample}.mmbwa.bam'
    """
}