import polars as pl

metadata = pl.read_csv(config['input_files']["metadata"], separator="\t")

samples = metadata["sample"].to_list()
read1 = metadata["R1"].to_list()
FIRST_SAMPLE = read1[0]
read2 = metadata["R2"].to_list()

READ_LOOKUP = dict(zip(samples, zip(read1, read2)))

STRAND_MAP = {
    "reverse": "--rna-strandness RF",
    "forward": "--rna-strandness FR",
    "unstranded": "--rna-strandness unstranded"
}

