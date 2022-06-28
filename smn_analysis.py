#!/usr/env python3

"""Counts and prints the number of reads that overlap exon 7 of SMN1 and SMN2 in the input CRAM or BAM file. 
"""
import os
import pysam
import sys
import pandas as pd

# parse arguments
reference_fasta_path = sys.argv[1]
cram_paths = sys.argv[2:]
output_file_location = '/Users/rakshya/SMA.2022/smn_analysis/count.tsv'


output_records = []
for cram_path in cram_paths:
    cram = pysam.AlignmentFile(cram_path, 'rc', reference_filename=reference_fasta_path)
    SMN1_exon7_read_count = cram.count(region="chr5:70924941-70953012")
    SMN2_exon7_read_count = cram.count(region="chr5:70049523-70077595")

    cram_filename = os.path.basename(cram_path)
    cram_filename = cram_filename.replace(".cram", "")

    record = {
        'SampleID': cram_filename,
        'SMN1': SMN1_exon7_read_count,
        'SMN2': SMN2_exon7_read_count,
        'SMN1 + SMN2': SMN1_exon7_read_count + SMN2_exon7_read_count,
    }

    output_records.append(record)

    cram.close()


df = pd.DataFrame(output_records)
df.to_csv(output_file_location, sep='\t', index=False)

print("No of cram files: ", len(output_records), " output in", output_file_location)

