from pathlib import Path

rule rsem_index:
    output:
        signal="results/rsem_index/.continue",
        index_dir=directory("results/rsem_index"),
    params:
        fasta=config['input_files']['reference_fasta'],
        gtf=config['input_files']['gtf_file'],
        transcript_to_gene_map=(f"--transcript-to-gene-map {Path(config['input_files']['knownIsoforms']).name}", config['input_files']['knownIsoforms']) if config['input_files'].get('knownIsoforms') else "",
        tmpdir=config['tmpdir'],
    container: config['container']
    log: "logs/rsem_index.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Creating rsem index"

        OUTDIR={output.index_dir}

        TEMP_FASTA={params.tmpdir}/$(uuidgen).fa
        TEMP_GTF={params.tmpdir}/$(uuidgen).gtf
        cp {params.transcript_to_gene_map[1]} {params.tmpdir}
        cp {params.fasta} $TEMP_FASTA
        cp {params.gtf} $TEMP_GTF

        cd {params.tmpdir}
        mkdir -p $OUTDIR

        rsem-prepare-reference \
            --gtf $TEMP_GTF \
            {params.transcript_to_gene_map[0]} \
            -p {resources.threads} \
            $TEMP_FASTA \
            {output.index_dir}/rsem

        rm -rf $TEMP_FASTA
        rm -rf $TEMP_GTF
        rm -rf $(basename {params.transcript_to_gene_map[1]})

        cd -

        mv {params.tmpdir}/$OUTDIR/* $OUTDIR

        echo "$(date): Finished creating rsem index"
        touch {output.signal}
        """

rule rsem_quant:
    input:
        signal=rules.rsem_index.output.signal,
        index_dir=rules.rsem_index.output.index_dir,
        anno_bam="results/align/{sample}/{sample}_anno.bam",
    output:
        outdir=directory("results/rsem_quant/{sample}"),
        signal="results/rsem_quant/{sample}/.continue",
    params:
        sample=lambda wildcards: wildcards.sample,
        endedness=config['sequencing']['endedness'],
        strand_direction=config['sequencing']['strand_direction'],
        rsem_seed=config['sequencing']['rsem_seed'],
        tmpdir=config['tmpdir'],
    container: config['container']
    log: "logs/rsem_quant/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        ABS_LOG_PATH=$(realpath {log})
        echo "$(date): Running rsem quantification"

        TARGET_DIR={params.tmpdir}/{output.outdir}

        TEMP_INDEX={params.tmpdir}/$(uuidgen)
        mkdir -p $TEMP_INDEX
        cp {input.anno_bam} {params.tmpdir}
        cp {input.index_dir}/* $TEMP_INDEX

        mkdir -p $TARGET_DIR
        cd $TARGET_DIR
    
        rsem.py \
            --rsem_index $TEMP_INDEX/rsem \
            --anno_bam {params.tmpdir}/$(basename {input.anno_bam}) \
            --logfile $ABS_LOG_PATH \
            --endedness {params.endedness} \
            --read_strand {params.strand_direction} \
            --rnd_seed {params.rsem_seed} \
            --ncpus {resources.threads} \
            --ramGB {resources.mem_gb}
        
        rm -rf $TEMP_INDEX
        rm -rf {params.tmpdir}/$(basename {input.anno_bam})
        
        cd -

        mv {params.tmpdir}/{output.outdir}/* {output.outdir}

        echo "$(date): Finished rsem quantification"
        touch {output.signal}
        """