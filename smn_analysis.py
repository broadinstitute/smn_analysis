#!/usr/env python

""" Read CRAM files"""
""" The cram files contains the data of SMN mutations diagnosed samples. """
import pysam
HUMAN_G = '/Users/rakshya/SMA.2022/Homo_sapiens_assembly38.fasta'
CRAM_F = '/Users/rakshya/SMA.2022/MAAC057.cram'
cram = pysam.AlignmentFile(CRAM_F, 'rc', reference_filename=HUMAN_G)
SMN1 = cram.count(region="chr5:70924941-70953012")
print("SMN1 : ", SMN1)
SMN2 = cram.count(region="chr5:70049523-70077595")
print("SMN2 : ", SMN2)
cram.close()

