rule rsem_index:
    output:
        rsem_index=directory("results/rsem_index"),
    params:
        fasta=config['input_files']['reference_fasta'],
        gtf=config['input_files']['gtf_file'],
        transcript_to_gene_map=f"--transcript-to-gene-map {config['input_files']['knownIsoforms']}" if config['input_files'].get('knownIsoforms') else "",
    container: config['container']
    log: "logs/rsem_index.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Creating rsem index"
        mkdir -p {output.rsem_index}
    
        rsem-prepare-reference \
            --gtf {params.gtf} \
            {params.transcript_to_gene_map} \
            --star \
            -p {resources.threads} \
            {params.fasta} \
            {output.rsem_index}/rsem

        echo "$(date): Finished creating rsem index"
        
        """

        # --transcript-to-gene-map {params.transcript_to_gene_map} \

rule rsem_quant:
    input:
        rsem_index=rules.rsem_index.output.rsem_index,
        anno_bam="results/align/{sample}/{sample}_anno.bam",
    output:
        outdir=directory("results/rsem_quant/{sample}"),
        signal="results/rsem_quant/{sample}/.continue",
    params:
        sample=lambda wildcards: wildcards.sample,
        endedness=config['sequencing']['endedness'],
        strand_direction=config['sequencing']['strand_direction'],
        rsem_seed=config['sequencing']['rsem_seed'],
    container: config['container']
    log: "logs/rsem_quant/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Running rsem quantification"

        SNAKE_WORKDIR=$(pwd)
    
        cd {output.outdir}
        $SNAKE_WORKDIR/workflow/rules/scripts/rsem.py \
            --rsem_index $SNAKE_WORKDIR/{input.rsem_index}/rsem \
            --anno_bam $SNAKE_WORKDIR/{input.anno_bam} \
            --logfile $SNAKE_WORKDIR/{log} \
            --endedness {params.endedness} \
            --read_strand {params.strand_direction} \
            --rnd_seed {params.rsem_seed} \
            --ncpus {resources.threads} \
            --ramGB {resources.mem_gb}
        cd -

        echo "$(date): Finished rsem quantification"
        touch {output.signal}
        """