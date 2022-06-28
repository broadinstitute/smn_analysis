#!/usr/env python3

"""Counts and prints the number of reads that overlap exon 7 of SMN1 and SMN2 in the input CRAM or BAM file. 
"""
import pysam
import sys
import pandas as pd

cram = pysam.AlignmentFile(sys.argv[2], 'rc', reference_filename=sys.argv[1])
SMN1 = cram.count(region="chr5:70924941-70953012")
SMN2 = cram.count(region="chr5:70049523-70077595")
file_location = '/Users/rakshya/SMA.2022/smn_analysis/count.tsv'
#parsing the cram filename
cram_filename = (sys.argv[2].split('/')[-1]).split('.')[0]
headerlist = ['SampleID','SMN1','SMN2','SMN1 + SMN2']
#inserting the data in count.tsv file
df = pd.DataFrame({'SampleID': [cram_filename], 'SMN1': [str(SMN1)], 'SMN2': [str(SMN2)], 'SMN1 + SMN2': [str(SMN1 + SMN2)]})
df.to_csv(file_location, mode = 'a', sep='\t', index=False)
cram.close()

"""//Output// 
SampleID        SMN1    SMN2    SMN1 + SMN2
MAAC057 7088    7438    14526
SampleID        SMN1    SMN2    SMN1 + SMN2
MANT033 7642    9056    16698
"""