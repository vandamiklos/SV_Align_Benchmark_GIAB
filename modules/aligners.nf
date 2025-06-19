process AlignReads {

    tag "$aligner"

    input:
    val aligner
    path ref
    path reads
    path sq_lines
    val threads
    path outdir

    output:
    path "${outdir}/${aligner}.*"

    script:
    def minimap2Preset = (aligner == 'minimap2') ? (params.platform == 'ont' ? "map-ont" : "map-hifi") : ""

    def aln_cmd = [
        'bwa': """(time bwa mem -t $threads $ref $reads) 2>> $outdir/bwa.log | \
                  samtools view -bh - | samtools sort -o $outdir/bwa.bam && samtools index $outdir/bwa.bam""",

        'minimap2': """(time minimap2 -ax $minimap2Preset -t $threads $ref $reads) 2>> $outdir/minimap2.log | \
                       samtools view -bh - | samtools sort -o $outdir/minimap2.bam && samtools index $outdir/minimap2.bam""",

        'last': """head -n 100 $reads > $outdir/training_set.fq && \
                   last-train -Q0 $ref $outdir/training_set.fq > $outdir/trained_parameters.txt && \
                   lastal --split -p $outdir/trained_parameters.txt -P $threads -Q0 $ref $reads | \
                   maf-convert sam > $outdir/lastalsplit.sam && \
                   samtools view -H $outdir/lastalsplit.sam > $outdir/old_header.sam && \
                   cat $outdir/old_header.sam $sq_lines > $outdir/new_header.sam && \
                   cat $outdir/new_header.sam $outdir/lastalsplit.sam | \
                   samtools view -bh - | samtools sort -o $outdir/lastalsplit.bam && \
                   samtools index $outdir/lastalsplit.bam""",

        'ngmlr': """ngmlr -t $threads -r $ref -q $reads | \
                    samtools view -bh - | samtools sort -o $outdir/ngmlr.bam && samtools index $outdir/ngmlr.bam""",

        'vacmap': """vacmap -ref $ref -read $reads -mode S -t $threads | \
                     samtools sort -o $outdir/vacmap.bam && samtools index $outdir/vacmap.bam""",

        'mmbwa': """echo 'MMBWA not implemented yet' > $outdir/mmbwa.log"""
    ][aligner]

    """
    mkdir -p $outdir
    $aln_cmd
    """
}
