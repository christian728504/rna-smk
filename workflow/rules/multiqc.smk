rule multiqc:
    input:
        expand("results/rsem_quant/{sample}.continue", sample=samples),
    output:
        datadir=directory("results/multiqc/multiqc_data"),
        html_report="results/multiqc/multiqc_report.html",
        signal="results/multiqc/.continue",
    params: 
        multiqc_conf=f"--multiqc-config {config["input_files"]['multiqc_config']}" if config['input_files'].get('multiqc_config') else "",
        outdir="results"
    container: config['container']
    log: f"logs/mutliqc.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date --iso=minutes): Started multiqc"

        python workflow/rules/scripts/run_multiqc.py {params.multiqc_conf} --outdir {params.outdir} --overwrite

        echo "$(date --iso=minutes): Finished multiqc"
        touch {output.signal}
        """
