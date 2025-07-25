def get_bam_to_signals_rule():
    if config['sequencing']['is_stranded']:
        return {
            'output': {
                'minus_all': "results/signals/{sample}/{sample}_genome_minusAll.bw",
                'minus_unique': "results/signals/{sample}/{sample}_genome_minusUniq.bw",
                'plus_all': "results/signals/{sample}/{sample}_genome_plusAll.bw",
                'plus_unique': "results/signals/{sample}/{sample}_genome_plusUniq.bw",
            }
        }
    else:
        return {
            'output': {
                'unstranded_all': "results/signals/{sample}/{sample}_genome_all.bw",
                'unstranded_unique': "results/signals/{sample}/{sample}_genome_uniq.bw",
            }
        }

rule bam_to_signals:
    input:
        chromsizes=rules.index.output.chromsizes,
        genome_bam="results/align/{sample}/{sample}_genome.bam"
    output:
        **get_bam_to_signals_rule()['output'],
        outdir=directory("results/signals/{sample}"),
        signal="results/signals/{sample}/.continue",
    params:
        sample=lambda wildcards: wildcards.sample,
        is_stranded=config['sequencing']['is_stranded'],
        tmpdir=config['tmpdir'],
    container: config['container']
    log: "logs/signals/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        ABS_LOG_PATH=$(realpath {log})
        echo "$(date): Running bam to singal conversion"

        OUTDIR={output.outdir}

        cp {input.genome_bam} {params.tmpdir}
        TEMP_CHROMSIZES={params.tmpdir}/$(uuidgen).txt
        cp {input.chromsizes} $TEMP_CHROMSIZES

        cd {params.tmpdir}
        mkdir -p $OUTDIR
        
        samtools sort -@ {resources.threads} -m 1G -o {output.outdir}/temp_sorted.bam $(basename {input.genome_bam})

        bam_to_signals.py \
            --logfile $ABS_LOG_PATH \
            --threads {resources.threads} \
            --output_dir {output.outdir} \
            --bamfile {output.outdir}/temp_sorted.bam \
            --chrom_sizes $TEMP_CHROMSIZES \
            --is_stranded {params.is_stranded} \
            --bamroot {params.sample}_genome

        rm -rf {output.outdir}/temp_sorted.bam
        rm -rf $TEMP_CHROMSIZES
        rm -rf $(basename {input.genome_bam})

        cd -

        mv {params.tmpdir}/$OUTDIR/* $OUTDIR

        echo "$(date): Finished bam to singal conversion"
        touch {output.signal}
        """