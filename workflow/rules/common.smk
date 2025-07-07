import polars as pl
import random

metadata = pl.read_csv(config['input_files']["metadata"], separator="\t")

samples = metadata["sample"].to_list()
read1 = metadata["R1"].to_list()
RANDOM_SAMPLE = random.choice(read1)
read2 = metadata["R2"].to_list()

READ_LOOKUP = dict(zip(samples, zip(read1, read2)))