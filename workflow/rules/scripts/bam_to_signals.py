#!/usr/bin/env python3
"""
Script to run bam to signals (bigwigs) step in ENCODE rna-seq-pipeline
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import logging
import shlex
import subprocess
import sys
import os

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)

STAR_COMMAND = """STAR --runMode inputAlignmentsFromBAM \
                --inputBAMfile {input_bam} \
                --outWigType bedGraph \
                --outWigStrand {strandedness} \
                --outWigReferencesPrefix chr
                --outFileNamePrefix {output_dir}/"""

def main(args):
    print(args)
    star_return_code = call_star(args.bamfile, args.is_stranded, args.output_dir)

    try:
        assert star_return_code == 0
    except AssertionError:
        logger.exception("Building bedGraph had a problem, most likely out of memory.")
        sys.exit(1)
        
    kwargs = {
            "chrom_sizes": args.chrom_sizes,
            "threads": args.threads,
        }
    
    if args.is_stranded:
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.UniqueMultiple.str1.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_minusAll.bw"),
            **kwargs,
        )
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.Unique.str1.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_minusUniq.bw"),
            **kwargs,
        )
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.UniqueMultiple.str2.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_plusAll.bw"),
            **kwargs,
        )
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.Unique.str2.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_plusUniq.bw"),
            **kwargs,
        )
    else:
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.UniqueMultiple.str1.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_all.bw"),
            **kwargs,
        )
        call_bg_to_bw(
            os.path.join(args.output_dir, "Signal.Unique.str1.out.bg"),
            os.path.join(args.output_dir, args.bamroot + "_uniq.bw"),
            **kwargs,
        )


def call_star(input_bam, is_stranded, output_dir):
    
    if is_stranded:
        strandedness = "stranded"
    else:
        strandedness = "unstranded"
    
    command = STAR_COMMAND.format(
        input_bam=input_bam, strandedness=strandedness.capitalize(), output_dir=output_dir
    )
    logger.info("Running STAR command %s", command)
    return_code = subprocess.call(shlex.split(command))
    return return_code


def call_bg_to_bw(input_bg, out_fn, chrom_sizes, threads):
    # sort bedgraph
    sorted_bg = input_bg.replace('.bg', '.sorted.bg')
    
    bedgraph_cmd = "bedtools sort -i {input_bg}".format(
        input_bg=input_bg
    )
    
    logger.info("Sorting bedgraph: %s", bedgraph_cmd)
    
    with open(sorted_bg, "w") as f:
        subprocess.call(shlex.split(bedgraph_cmd), stdout=f)
    # make bigwig
    command = "bigtools bedgraphtobigwig --nthreads {nthreads} {sorted_bg} {chrom_sizes} {out_fn}".format(
        sorted_bg=sorted_bg, chrom_sizes=chrom_sizes, out_fn=out_fn, nthreads=threads
    )
    
    logger.info("Building bigWig: %s", command)
    
    subprocess.call(shlex.split(command))
    os.remove(input_bg)
    os.remove(sorted_bg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--logfile", help="Path to log file.", required=True)
    parser.add_argument("--output_dir", type=str, help="Path to output directory.", required=True)
    parser.add_argument("--threads", type=int, help="Number of threads for `bigtools bedgraphtobigwig`", default=1) 
    parser.add_argument("--bamfile", type=str, help="Input bam")
    parser.add_argument(
        "--chrom_sizes", type=str, help="chromosome sizes file the input bam"
    )
    parser.add_argument("--is_stranded", type=bool)
    parser.add_argument(
        "--bamroot",
        type=str,
        help="""
             Root name for output bams. For example out_bam
             will create out_bam_genome.bam and out_bam_anno.bam
             """,
        default="out_bam",
    )
    args = parser.parse_args()
    
    filehandler = logging.FileHandler(args.logfile)
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(formatter)
    logger.addHandler(filehandler)
    
    main(args)
