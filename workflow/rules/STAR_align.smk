rule index:
    output:
        chromsizes="results/index/chrNameLength.txt",
        signal="results/index/.continue",
        index_dir=directory("results/index"),
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
        echo "$(date): Creating STAR index"

        OUTDIR={output.index_dir}

        TEMP_FASTA={params.tmpdir}/$(uuidgen).fa
        TEMP_GTF={params.tmpdir}/$(uuidgen).gtf
        cp {params.sample} {params.tmpdir}
        cp {params.fasta} $TEMP_FASTA
        cp {params.gtf} $TEMP_GTF

        cd {params.tmpdir}
        mkdir -p $OUTDIR

        AVG_READ_LENGTH=$(seqkit stats $(basename {params.sample}) | awk 'NR>1 {{print int($6)}}')

        STAR \
            --runThreadN {resources.threads} \
            --runMode genomeGenerate \
            --genomeDir {output.index_dir} \
            --genomeFastaFiles $TEMP_FASTA \
            --sjdbGTFfile $TEMP_GTF \
            --sjdbOverhang $(($AVG_READ_LENGTH - 1))

        rm -rf $TEMP_FASTA
        rm -rf $TEMP_GTF
        rm -rf $(basename {params.sample})

        cd -

        mv {params.tmpdir}/$OUTDIR/* $OUTDIR

        echo "$(date): Finished creating STAR index"
        touch {output.signal}
        """

rule align:
    input:
        signal=rules.index.output.signal,
        index_dir=rules.index.output.index_dir,
        r1="results/processed/{sample}/{sample}_R1.fastq.gz",
        r2="results/processed/{sample}/{sample}_R2.fastq.gz"
    output:
        genome_bam="results/align/{sample}/{sample}_genome.bam",
        anno_bam="results/align/{sample}/{sample}_anno.bam",
        genome_flagstat="results/align/{sample}/{sample}_genome_flagstat.txt",
        anno_flagstat="results/align/{sample}/{sample}_anno_flagstat.txt",
        log="results/align/{sample}/{sample}_Log.final.out",
        genome_flagstat_json="results/align/{sample}/{sample}_genome_flagstat.json",
        anno_flagstat_json="results/align/{sample}/{sample}_anno_flagstat.json",
        log_json="results/align/{sample}/{sample}_Log.final.json",
        outdir=directory("results/align/{sample}"),
        signal="results/align/{sample}/.continue"
    params:
        sample=lambda wildcards: wildcards.sample,
        endedness=config['sequencing']['endedness'],
        tmpdir=config['tmpdir'],
    container: config['container']
    log: "logs/align/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        ulimit -n 65636
        ABS_LOG_PATH=$(realpath {log})
        ABS_SCRIPT_PATH=$(realpath workflow/rules/scripts/align.py)
        echo "$(date): Running STAR alignment"

        OUTDIR={output.outdir}

        cp {input.r1} {input.r2} {params.tmpdir}
        TEMP_INDEX={params.tmpdir}/$(uuidgen)
        mkdir -p $TEMP_INDEX
        cp {input.index_dir}/* $TEMP_INDEX

        cd {params.tmpdir}
        mkdir -p $OUTDIR

        $ABS_SCRIPT_PATH \
            --logfile $ABS_LOG_PATH \
            --output_dir {output.outdir} \
            --fastqs_R1 $(basename {input.r1}) \
            --fastqs_R2 $(basename {input.r2}) \
            --bamroot {params.sample} \
            --indexdir $TEMP_INDEX \
            --endedness {params.endedness} \
            --ramGB {resources.mem_gb} \
            --ncpus {resources.threads} 
        
        samtools quickcheck {output.genome_bam} {output.anno_bam} && \
            echo "$(date): samtools quickcheck OK" || echo "$(date): samtools quickcheck FAIL"

        rm -rf $(basename {input.r1})
        rm -rf $(basename {input.r2})
        rm -rf $TEMP_INDEX

        cd -

        mv {params.tmpdir}/$OUTDIR/* $OUTDIR

        echo "$(date): Finished STAR alignment"
        touch {output.signal}
        """
