rule star_index:
    output:
        chromsizes="results/genome/chrNameLength.txt",
        genome_dir=directory("results/genome")
    params:
        sample=RANDOM_SAMPLE,
        fasta=config['input_files']['reference_fasta'],
        gtf=config['input_files']['gtf_file']
    container: config['container']
    log: "logs/index.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Creating STAR index"

        AVG_READ_LENGTH=$(seqkit stats {params.sample} | awk 'NR>1 {{print int($6)}}')

        STAR \
            --runThreadN {resources.threads} \
            --runMode genomeGenerate \
            --genomeDir {output.genome_dir} \
            --genomeFastaFiles {params.fasta} \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang $(($AVG_READ_LENGTH - 1))

        echo "$(date): Finished creating STAR index"
        """

rule star_align:
    input:
        genome_dir=rules.star_index.output.genome_dir
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
        r1=lambda wildcards: READ_LOOKUP[wildcards.sample][0],
        r2=lambda wildcards: READ_LOOKUP[wildcards.sample][1],
        endedness=config['sequencing']['endedness']
    container: config['container']
    log: "logs/align/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Running STAR alignment"

        ./workflow/rules/scripts/align.py \
            --logfile {log} \
            --output_dir {output.outdir} \
            --fastqs_R1 {params.r1} \
            --fastqs_R2 {params.r2} \
            --bamroot {params.sample} \
            --indexdir {input.genome_dir} \
            --endedness {params.endedness} \
            --ramGB {resources.mem_gb} \
            --ncpus {resources.threads} 
        
        samtools quickcheck {output.genome_bam} {output.anno_bam} && \
            echo "$(date): samtools quickcheck OK" || echo "$(date): samtools quickcheck FAIL"

        echo "$(date): Finished STAR alignment"
        touch {output.signal}
        """