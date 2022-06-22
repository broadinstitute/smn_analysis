#!/usr/env python

""" Read CRAM files"""
import pysam
HUMAN_G = '/Users/rakshya/SMA.2022/Homo_sapiens_assembly38.fasta'
CRAM_FH = '/Users/rakshya/SMA.2022/MAAC057.cram'
cram = pysam.AlignmentFile(CRAM_FH, 'rc', reference_filename=None)
with pysam.FastxFile(HUMAN_G) as fh:
    for entry in fh:
        print(entry)
#SMN1
for read in cram.fetch(region="chr5:70924941-70953012", until_eof=True):
    print(read)
#SMN2
for read in cram.fetch(region="chr5:70049523-70077595", until_eof=True):
    print(read)
cram.close()

