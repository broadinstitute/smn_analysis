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
output_file_location = '/Users/rakshya/SMA.2022/smn_analysis/add-all-exon-count.tsv'
output_records = []
EXON_COORDINATES = [("chr5:70925087-70925184","chr5:69345439-69345593"),
                    ("chr5:70938839-70938910","chr5:69359242-69359313"),
                    ("chr5:70941389-70941508","chr5:69361792-69361911"),
                    ("chr5:70942358-70942558","chr5:69362761-69362961"),
                    ("chr5:70942718-70942870","chr5:69363121-69363273"),
                    ("chr5:70944658-70944753","chr5:69365062-69365157"),
                    ("chr5:70946066-70946176","chr5:69366468-69366578"),
                    ("chr5:70951941-70951994","chr5:69372348-69372401"),
                    ("chr5:70952439-70953015","chr5:69372846-69373419")]
for cram_path in cram_paths:
    cram = pysam.AlignmentFile(cram_path, 'rc', reference_filename=reference_fasta_path)
    #listing exon coordinates exon 1-8
    cram_filename = os.path.basename(cram_path)
    cram_filename = cram_filename.replace(".cram", "")
    #empty dictionary
    record = {}
    count = 0
    #updating the records dictionary
    record.update({"SAMPLEID": cram_filename})
    #looping to read exon coordinates from the list
    for SMN1_exon_region,SMN2_exon_region in EXON_COORDINATES:
        SMN1_exon_read_count = cram.count(region= SMN1_exon_region)
        SMN2_exon_read_count = cram.count(region= SMN2_exon_region)
        count+= 1
        exon = "SMN1 + SMN2 exon " + str(count) + " reads"
        exon_read = {exon : SMN1_exon_read_count + SMN2_exon_read_count}
        record.update(exon_read)

    output_records.append(record)
    cram.close()
df = pd.DataFrame(output_records)
df.to_csv(output_file_location, sep='\t', index=False)
print("No of cram files: ", len(output_records), " output in", output_file_location)

