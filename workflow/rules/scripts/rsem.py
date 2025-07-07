#!/usr/bin/env python3
"""
Script to run alignment(mapping) step in ENCODE rna-seq-pipeline
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import json
import logging
import os
import re
import shlex
import subprocess
import tarfile

import pandas as pd
from qc_utils import QCMetric, QCMetricRecord

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)

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


def strand_to_fwd_prob(strand):
    """Converts strand into a numeric value that RSEM understands.

    Args:
        strand: string 'forward', 'reverse', 'unstranded'

    Returns:
        numeric value corresponding the forward strand probability

    Raises:
        KeyError if strand is not 'forward', 'reverse' or 'unstranded'
    """
    conversion = {"forward": 1, "unstranded": 0.5, "reverse": 0}
    return conversion[strand]


def format_endedness(endedness):
    if endedness == "paired":
        return "--paired-end"
    else:
        return ""


def calculate_number_of_genes_detected(quant_tsv, threshold_of_detection=1):
    """Calculates number of rows where value on the column TPM is greater than
    the threshold_of_detection.

    Args:
        quant_tsv: filename of a .tsv of RSEM quants

    Returns:
        int number_of_genes_detected: number of entries > threshold_of_detection
        in TPM column
    """
    quants = pd.read_csv(quant_tsv, sep="\t")
    number_of_genes_detected = sum(quants["TPM"] > threshold_of_detection)
    return number_of_genes_detected


def main(args):
    remove_bam_from_end_re = re.compile("\.bam$")
    bam_root = remove_bam_from_end_re.sub("", os.path.basename(args.anno_bam))
    rsem_call = shlex.split(
        RSEM_COMMAND.format(
            rnd_seed=args.rnd_seed,
            ncpus=args.ncpus,
            ramGB=args.ramGB,
            fwd_prob=strand_to_fwd_prob(args.read_strand),
            paired_end=format_endedness(args.endedness),
            anno_bam=args.anno_bam,
            bam_root=bam_root,
            resem_index=args.rsem_index,
        )
    )
    logger.info("Running RSEM command %s", " ".join(rsem_call))
    subprocess.call(rsem_call)
    gene_quant_fn = str(bam_root) + "_rsem.genes.results"
    number_of_genes_detected = calculate_number_of_genes_detected(gene_quant_fn)
    number_of_genes_detected_dict = {
        "number_of_genes_detected": number_of_genes_detected
    }
    qc_record = QCMetricRecord()
    number_of_genes_QC = QCMetric(
        "number_of_genes_detected", number_of_genes_detected_dict
    )
    qc_record.add(number_of_genes_QC)

    with open(str(bam_root) + "_number_of_genes_detected.json", "w") as f:
        json.dump(qc_record.to_ordered_dict(), f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rsem_index", type=str, help="path to rsem index")
    parser.add_argument("--anno_bam", type=str, help="STAR alignment to annotation.")
    parser.add_argument("--logfile", type=str, help="path to log file")
    parser.add_argument("--endedness", type=str, choices=["paired", "single"])
    parser.add_argument(
        "--read_strand", type=str, choices=["forward", "reverse", "unstranded"]
    )
    parser.add_argument("--rnd_seed", type=int, help="random seed", default=42)
    parser.add_argument("--ncpus", type=int, help="number of cpus available")
    parser.add_argument("--ramGB", type=int, help="memory available in GB")
    args = parser.parse_args()
    
    filehandler = logging.FileHandler(args.logfile)
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(formatter)
    logger.addHandler(filehandler)
    
    main(args)
