#!/usr/env python

"""Counts and prints the number of reads that overlap exon 7 of SMN1 and SMN2 in the input CRAM or BAM file. 
"""
import pysam
REFERENCE_FASTA_PATH = '/Users/rakshya/SMA.2022/Homo_sapiens_assembly38.fasta'
CRAM_F = '/Users/rakshya/SMA.2022/MAAC057.cram'
cram = pysam.AlignmentFile(CRAM_F, 'rc', reference_filename=HUMAN_G)
SMN1 = cram.count(region="chr5:70924941-70953012")
print("SMN1 : ", SMN1)
#Output - SMN1 :  7088
SMN2 = cram.count(region="chr5:70049523-70077595")
print("SMN2 : ", SMN2)
#Output - SMN2 :  7438
cram.close()

