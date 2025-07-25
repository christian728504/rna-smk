---
theme:
  name: dark
  override:
    code:
      alignment: left
title: __A rework of the ENCODE DCC rna-seq-pipeline__
author: Christian S. Ramirez
---

RNA-seq workflow
===

First an we trim and QC the fastq files:

```bash
fastp -i $(basename {params.r1}) -I $(basename {params.r2}) \
              -o {output.r1} -O {output.r2} \
              -h {output.report} -j {output.json} \
              -w {resources.threads}
```

`fastp` is a will trim adapters and low quality bases from reads. It will also generate a report and a json file with key metrics about the reads, both before and after trimming. This will be visualized in the `multiqc` report.

<!-- end_slide -->

RNA-seq workflow
===

First an index is created for HISAT2:

```bash
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
```

<!-- end_slide -->

RNA-seq workflow
===

Also in STAR:

```bash
AVG_READ_LENGTH=$(seqkit stats $(basename {params.sample}) | awk 'NR>1 {{print int($6)}}')

STAR \
    --runThreadN {resources.threads} \
    --runMode genomeGenerate \
    --genomeDir {output.index_dir} \
    --genomeFastaFiles $TEMP_FASTA \
    --sjdbGTFfile $TEMP_GTF \
    --sjdbOverhang $(($AVG_READ_LENGTH - 1))
```

<!-- end_slide -->

RNA-seq workflow
===

Running the alignment in HISAT2:

```bash
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
```

<!-- end_slide -->

RNA-seq workflow
===

Continued:

```bash
amtools view -@ {resources.threads} -h -F 0x100 {output.genome_bam} > {output.outdir}/temp_primary.bam
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
```

<!-- end_slide -->

RNA-seq workflow
===

Also in STAR:

```bash
align.py \
    --logfile $ABS_LOG_PATH \
    --output_dir {output.outdir} \
    --fastqs_R1 $(basename {input.r1}) \
    --fastqs_R2 $(basename {input.r2}) \
    --bamroot {params.sample} \
    --indexdir $TEMP_INDEX \
    --endedness {params.endedness} \
    --ramGB {resources.mem_gb} \
    --ncpus {resources.threads} 
```

<!-- end_slide -->

RNA-seq workflow
===

Effective command in the `align.py` script:

```python
command_string = """STAR --genomeDir {indexdir} \
--readFilesIn {read1_fq_gz} {read2_fq_gz} \
--readFilesCommand zcat \
--runThreadN {ncpus} \
--genomeLoad NoSharedMemory \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMheaderCommentFile COfile.txt \
--outSAMheaderHD @HD VN:1.4 SO:coordinate \
--outSAMunmapped Within \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--sjdbScore 1 \
--limitBAMsortRAM {ramGB}000000000"""
```

<!-- end_slide -->

RNA-seq workflow
===

Then rsem for quantification. First we build the index:

```bash
 rsem-prepare-reference \
    --gtf $TEMP_GTF \
    {params.transcript_to_gene_map[0]} \
    -p {resources.threads} \
    $TEMP_FASTA \
    {output.index_dir}/rsem
```

<!-- end_slide -->

RNA-seq workflow
===

Running `rsem-calculate-expression`

```bash
rsem.py \
    --rsem_index $TEMP_INDEX/rsem \
    --anno_bam {params.tmpdir}/$(basename {input.anno_bam}) \
    --logfile $ABS_LOG_PATH \
    --endedness {params.endedness} \
    --read_strand {params.strand_direction} \
    --rnd_seed {params.rsem_seed} \
    --ncpus {resources.threads} \
    --ramGB {resources.mem_gb}
```

<!-- end_slide -->

RNA-seq workflow
===

Effective command in the `rsem.py` script:

```python
RSEM_COMMAND = """rsem-calculate-expression --bam \
--estimate-rspd \
--calc-ci \
--seed {rnd_seed} \
-p {ncpus} \
--no-bam-output \
--ci-memory {ramGB}000 \
--forward-prob {fwd_prob} \
{paired_end} \
{anno_bam} \
{resem_index} \
{bam_root}_rsem"""
```

<!-- end_slide -->

RNA-seq workflow
===

Effective command in the `rsem.py` script:

```python
RSEM_COMMAND = """rsem-calculate-expression --bam \
--estimate-rspd \
--calc-ci \
--seed {rnd_seed} \
-p {ncpus} \
--no-bam-output \
--ci-memory {ramGB}000 \
--forward-prob {fwd_prob} \
{paired_end} \
{anno_bam} \
{resem_index} \
{bam_root}_rsem"""
```

<!-- end_slide -->

RNA-seq workflow
===

Creating signal files:

For brevity, here are the effective commands from the `bam_to_signals.py` script:

```python
STAR_COMMAND = """STAR --runMode inputAlignmentsFromBAM \
                --inputBAMfile {input_bam} \
                --outWigType bedGraph \
                --outWigStrand {strandedness} \
                --outWigReferencesPrefix chr
                --outFileNamePrefix {output_dir}/"""
```

```python
def call_bg_to_bw(input_bg, out_fn, chrom_sizes, threads):
    sorted_bg = input_bg.replace('.bg', '.sorted.bg')
    bedgraph_cmd = "bedtools sort -i {input_bg}".format(input_bg=input_bg)
    subprocess.call(shlex.split(bedgraph_cmd), stdout=f)
    command = "bigtools bedgraphtobigwig --nthreads {nthreads} {sorted_bg} {chrom_sizes} {out_fn}".format(
        sorted_bg=sorted_bg, chrom_sizes=chrom_sizes, out_fn=out_fn, nthreads=threads
    )
    return(0)
```

<!-- end_slide -->

RNA-seq workflow
===

Creating signal files:

Lastly, we run `multiqc` to generate a nice report of all metrics:

```shell
multiqc {params.multiqc_conf} results --outdir results/multiqc
```

For the MOHD data, we have a multiqc config file to correcly identify sample wildcards in different files:

```yaml
extra_fn_clean_exts:
  - type: regex_keep
    pattern: ER\d+
```

<!-- end_slide -->

RNA-seq workflow
===

<!-- column_layout: [1, 3, 1] -->

<!-- column: 1 -->

Now to show the mermaid diagram as well as the multiqc report.

<!-- end_slide -->
