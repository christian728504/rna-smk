rule fastp:
    output:
        r1="results/processed/{sample}/{sample}_R1.fastq.gz",
        r2="results/processed/{sample}/{sample}_R2.fastq.gz",
        report="results/processed/{sample}/{sample}.html",
        json="results/processed/{sample}/{sample}.json",
        outdir=directory("results/processed/{sample}"),
    params:
        sample=lambda wildcards: wildcards.sample,
        r1=lambda wildcards: READ_LOOKUP[wildcards.sample][0],
        r2=lambda wildcards: READ_LOOKUP[wildcards.sample][1],
        tmpdir=config['tmpdir'],
    container: "docker://clarity001/wgs-smk:latest"
    log: "logs/processed/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Started fastp for sample: {params.sample}"
    
        OUTDIR=results/processed/{wildcards.sample}

        cp {params.r1} {params.r2} {params.tmpdir}

        cd {params.tmpdir}
        mkdir -p $OUTDIR

        fastp -i $(basename {params.r1}) -I $(basename {params.r2}) \
              -o {output.r1} -O {output.r2} \
              -h {output.report} -j {output.json} \
              -w {resources.threads}

        rm -rf $(basename {params.r1})
        rm -rf $(basename {params.r2})
        
        cd -

        mv {params.tmpdir}/$OUTDIR/* $OUTDIR
        
        echo "$(date): Completed fastp for sample: {params.sample}"
        """