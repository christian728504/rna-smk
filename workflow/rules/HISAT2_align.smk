rule index:
    output:
        chromsizes="results/index/chrNameLength.txt",
        signal="results/index/.continue",
        index_dir=directory("results/index"),
        splice_sites="results/index/genome.ss",
        exons="results/index/genome.exon",
    params:
        sample=FIRST_SAMPLE,
        fasta=config['input_files']['reference_fasta'],
        gtf=config['input_files']['gtf_file'],
        tmpdir=config['tmpdir'],
    container: config['container']
    log: "logs/index.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Creating HISAT2 index"

        OUTDIR={output.index_dir}

        TEMP_FASTA={params.tmpdir}/$(uuidgen).fa
        TEMP_GTF={params.tmpdir}/$(uuidgen).gtf
        cp {params.fasta} $TEMP_FASTA
        cp {params.gtf} $TEMP_GTF

        cd {params.tmpdir}
        mkdir -p $OUTDIR

        mudskipper index --dir-index {output.index_dir} --gtf $TEMP_GTF
        gtfextract exons --gtf $TEMP_GTF > {output.exons}
        gtfextract splice-sites --gtf $TEMP_GTF > {output.splice_sites}
        hisat2-build \
            -p {resources.threads} \
            --exon {output.exons} \
            --ss {output.splice_sites} \
            $TEMP_FASTA \
            {output.index_dir}/hisat2 > /dev/null
        samtools faidx -@ {resources.threads} $TEMP_FASTA -o - | cut -f 1,2 > {output.chromsizes}

        rm -rf $TEMP_FASTA
        rm -rf $TEMP_GTF

        cd -

        mv {params.tmpdir}/$OUTDIR/* $OUTDIR

        echo "$(date): Finished creating HISAT2 index"
        touch {output.signal}
        """

rule align:
    input:
        signal=rules.index.output.signal,
        index_dir=rules.index.output.index_dir,
        r1="results/processed/{sample}/{sample}_R1.fastq.gz",
        r2="results/processed/{sample}/{sample}_R2.fastq.gz",
    output:
        genome_bam="results/align/{sample}/{sample}_genome.bam",
        alignment_summary="results/align/{sample}/{sample}_alignment_summary.txt",
        anno_bam="results/align/{sample}/{sample}_anno.bam",
        genome_flagstat_json="results/align/{sample}/{sample}_genome_flagstat.json",
        anno_flagstat_json="results/align/{sample}/{sample}_anno_flagstat.json",
        outdir=directory("results/align/{sample}"),
        signal="results/align/{sample}/.continue"
    params:
        sample=lambda wildcards: wildcards.sample,
        endedness=config['sequencing']['endedness'],
        strand_direction=STRAND_MAP.get(config['sequencing']['strand_direction'], ""),
        tmpdir=config['tmpdir'],
    container: config['container']
    log: "logs/align/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Running HISAT2 alignment"

        OUTDIR={output.outdir}

        cp {input.r1} {input.r2} {params.tmpdir}
        TEMP_INDEX={params.tmpdir}/$(uuidgen)
        mkdir -p $TEMP_INDEX
        cp {input.index_dir}/* $TEMP_INDEX

        cd {params.tmpdir}
        mkdir -p $OUTDIR

        hisat2 \
            --threads {resources.threads} \
            --rg-id={params.sample} \
            {params.strand_direction} \
            -k 20 \
            --min-intronlen 20 \
            --max-intronlen 1000000 \
            --no-mixed \
            --no-discordant \
            --no-softclip \
            -x $TEMP_INDEX/hisat2 \
            --summary-file {output.alignment_summary} \
            -1 $(basename {input.r1}) \
            -2 $(basename {input.r2}) \
            -S {output.outdir}/temp_genome.sam
        samtools view -@ {resources.threads} -b -S {output.outdir}/temp_genome.sam > {output.genome_bam}

        samtools view -@ {resources.threads} -h -F 0x100 {output.genome_bam} > {output.outdir}/temp_primary.bam
        mudskipper bulk \
            --threads {resources.threads} \
            --max_mem_mb {resources.mem_mb} \
            --index $TEMP_INDEX \
            --alignment {output.outdir}/temp_primary.bam \
            --out {output.outdir}/temp_anno.bam \
            --skip \
            --shuffle > /dev/null
        filter-indels.py \
            {output.outdir}/temp_anno.bam \
            {output.outdir}/temp_noindels.bam \
            --threads {resources.threads}

        convert-sam-for-rsem \
            --num-threads {resources.threads} \
            --memory-per-thread $(( ({resources.mem_mb} / {resources.threads}) * 80 / 100 ))M \
            {output.outdir}/temp_noindels.bam \
            {output.outdir}/{wildcards.sample}_anno

        samtools flagstat \
            -@ {resources.threads} \
            --output-fmt json {output.genome_bam} > {output.genome_flagstat_json}
        samtools flagstat \
            -@ {resources.threads} \
            --output-fmt json {output.anno_bam} > {output.anno_flagstat_json}
        samtools quickcheck {output.genome_bam} {output.anno_bam} && \
            echo "$(date): samtools quickcheck OK" || echo "$(date): samtools quickcheck FAIL"
        
        rm -rf {output.outdir}/temp_* {output.outdir}/*tmp*
        rm -rf $(basename {input.r1})
        rm -rf $(basename {input.r2})
        rm -rf $TEMP_INDEX

        cd -

        mv {params.tmpdir}/$OUTDIR/* $OUTDIR

        echo "$(date): Finished HISAT2 alignment"
        touch {output.signal}
        """
