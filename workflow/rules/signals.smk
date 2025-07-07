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
        chromsizes=rules.star_index.output.chromsizes,
        genome_bam="results/align/{sample}/{sample}_genome.bam"
    output:
        **get_bam_to_signals_rule()['output'],
        outdir=directory("results/signals/{sample}"),
        signal="results/signals/{sample}/.continue",
    params:
        sample=lambda wildcards: wildcards.sample,
        is_stranded=config['sequencing']['is_stranded'],
    container: config['container']
    log: "logs/signals/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Running STAR bam to singal conversion"
        
        ./workflow/rules/scripts/bam_to_signals.py \
            --logfile {log} \
            --threads {resources.threads} \
            --output_dir {output.outdir} \
            --bamfile {input.genome_bam} \
            --chrom_sizes {input.chromsizes} \
            --is_stranded {params.is_stranded} \
            --bamroot {params.sample}_genome

        echo "$(date): Finished STAR bam to singal conversion"
        touch {output.signal}
        """