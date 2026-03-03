#!/usr/bin/env python

import argparse
from collections import defaultdict
import subprocess
import pysam
import tempfile
import sys
import os
import shutil

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def main(args):
    if args.src:
        try:
            editor = os.environ.get("EDITOR")
            if not editor:
                editor = "vi"
            subprocess.check_call([editor, sys.argv[0]])
            sys.exit(0)
        except subprocess.CalledProcessError:
            print("Error: Could not open file in editor.")
    
    in_path = args.input
    out_path = args.output
    threads = args.threads
    
    print("Filtering out indels from %s" % in_path)
    
    in_bam = pysam.AlignmentFile(in_path, 'rb')
    if not in_bam.has_index():
        print("No index found for input bam. Preparing index...")
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            temp_path = tmp.name
            shutil.copy(in_path, temp_path)
            
            print("Sorting input bam...")
            subprocess.check_call(["samtools", "sort", "-@", str(threads), "-o", in_path, temp_path])
        
        os.unlink(temp_path)
        print("Indexing input bam...")
        subprocess.check_call(["samtools", "index", "-@", str(threads), in_path])
        in_bam.close()
    
    print("Writing output bam...")
    in_bam = pysam.AlignmentFile(in_path, 'rb')
    out_bam = pysam.AlignmentFile(out_path, 'wb', template=in_bam)
    for read1, read2 in read_pair_generator(in_bam):
        if "I" in read1.cigarstring + read2.cigarstring or "D" in read1.cigarstring + read2.cigarstring:
            continue
        out_bam.write(read1)
        out_bam.write(read2)
    out_bam.close()
    in_bam.close()

if __name__ == "__main__":
    os.path.abspath(sys.argv[0])
    parser = argparse.ArgumentParser(
        prog="filter-indels",
        description="Tool for filtering out alignment records (in pairs) containing insertions (I) or deletions (D).",
    )
    parser.add_argument(
        "-t", "--threads", help="Number of threads provided to samtools", metavar="<THREADS>", type=int, default=1 
    )
    parser.add_argument(
        "--src", help="View source code", action="store_true"
    )
    parser.add_argument(
        "input", help="Path to input bam. Pass '-' to read from stdin.", nargs="?", metavar="<IN.BAM>", type=str,
    )
    parser.add_argument(
        "output", help="Path to output bam. Pass '-' to write to stdout.", nargs="?", metavar="<OUT.BAM>", type=str,
    )
    args = parser.parse_args()
    main(args)