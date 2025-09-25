#!/usr/bin/env python3

import sys
from Bio import SeqIO
import pandas as pd

reference_file = sys.argv[1]
bed_file = sys.argv[2]
patched_bed_file_name = f"{bed_file.replace('.bed', '.patched')}.bed"

print(f"Reading input files REF: '{reference_file}', BED: '{bed_file}' ...")

ref_seq = next(SeqIO.parse(reference_file, "fasta"))
bed_df = pd.read_csv(bed_file, sep="\t", names = ['chrom', 'chromStart', 'chromEnd', 'primer-name', 'pool', 'strand'])

if len(bed_df.columns) > 6:
    print('BED file is not in deprecated format and does not have to be patched for artic, how did you end up in this script?')
    quit()

print('Patching BED file primer sequences ...')

bed_df['primer-sequence'] = bed_df.apply(
    lambda row: str(ref_seq.seq[row['chromStart']:row['chromEnd']]) if row['strand'] == "+" else str(ref_seq.seq[row['chromStart']:row['chromEnd']].reverse_complement()),
    axis=1
)

print(f"Writing patched BED file '{patched_bed_file_name}' ...")

bed_df.to_csv(patched_bed_file_name, sep = '\t', header = False, index = False)

print(f"Done.")
