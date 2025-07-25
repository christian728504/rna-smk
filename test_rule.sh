#!/bin/bash


exec >> logs/align/ER100002.log 2>&1
ABS_LOG_PATH=$(realpath logs/align/ER100002.log)
ABS_SCRIPT_PATH=$(realpath ./workflow/rules/scripts/align.py)
echo "$(date): Running STAR alignment"

OUTDIR=results/align/ER100002

cp results/processed/ER100002/ER100002_R1.fastq.gz results/processed/ER100002/ER100002_R2.fastq.gz /tmp
TEMP_INDEX=/tmp/$(uuidgen)
mkdir -p $TEMP_INDEX
cp results/index/* $TEMP_INDEX

cd /tmp
mkdir -p $OUTDIR

$ABS_SCRIPT_PATH             --logfile $ABS_LOG_PATH             --output_dir results/align/ER100002             --fastqs_R1 $(basename results/processed/ER100002/ER100002_R1.fastq.gz)             --fastqs_R2 $(basename results/processed/ER100002/ER100002_R2.fastq.gz)             --bamroot ER100002             --indexdir $TEMP_INDEX             --endedness paired             --ramGB 125 --ncpus 22

samtools quickcheck results/align/ER100002/ER100002_genome.bam results/align/ER100002/ER100002_anno.bam &&             echo "$(date): samtools quickcheck OK" || echo "$(date): samtools quickcheck FAIL"

# rm -rf $(basename results/processed/ER100001/ER100002_R1.fastq.gz)
# rm -rf $(basename results/processed/ER100002/ER100002_R2.fastq.gz)
# rm -rf $TEMP_INDEX

cd -

# mv /tmp/$OUTDIR/* $OUTDIR

echo "$(date): Finished STAR alignment"
# touch results/align/ER100002/.continue
