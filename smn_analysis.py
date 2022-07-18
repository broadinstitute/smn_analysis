#!/usr/env python3

"""This script counts the number of reads that overlap exons of SMN1 and SMN2 in the input CRAM file."""

import argparse
import os
import pysam
import pandas as pd
import matplotlib.pyplot as plt
SMN1_C840_POSITION_1BASED = 70951946
SMN2_C840_POSITION_1BASED = 70076526

def parse_args():
    """Parse and return command-line args"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file path")
    parser.add_argument("-o", "--output-tsv", default="smn-exon-counts.tsv", help="Output tsv file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Whether to print extra stats and output")
    parser.add_argument("cram_paths", nargs="+", help="One or more CRAM (or BAM) file paths")
    args = parser.parse_args()

    # make sure the input files exist
    for path in [args.reference_fasta] + args.cram_paths:
        if not os.path.isfile(path):
            parser.error(f"File not found: {path}")

    return args


def compute_total_reads(cram_path):
    """Compute the total number of reads in the given cram file using samtools idxstat.
    Args:
        cram_path (str): local file path
    Return:
        int: The total number of reads
    """
    cram_filename = os.path.basename(cram_path)
    idxstat_filename = f"{cram_filename}.idxstat.txt"
    os.system(f"samtools idxstat {cram_path}> {idxstat_filename}")
    idxstat_df = pd.read_table(idxstat_filename, names=[
        "chrom", "chrom_length", "mapped_read_counts", "unmapped_read_counts",
    ])
    os.remove(idxstat_filename)

    total_reads = sum(idxstat_df.mapped_read_counts) + sum(idxstat_df.unmapped_read_counts)

    return total_reads

def nucleotide_count(cram,start,end):
    smn_nucleotide = {"T":0,"C":0,"G":0,"A":0,"N":0}
    for pileup_column in cram.pileup(contig="chr5", start=start, end=end, truncate=True):
        for pileupread in pileup_column.pileups:
            val = pileupread.alignment.query_sequence[pileupread.query_position]
            val = val.upper()
            if val in smn_nucleotide:
                smn_nucleotide[val] += 1
            else:
                raise Exception(f'Unexpected base "{val}" found in {pileupread}')
    return smn_nucleotide

def main():
    args = parse_args()

    if args.verbose:
        print(f"Processing {len(args.cram_paths)} CRAM files")

    output_records = []
    for cram_path in args.cram_paths:
        cram = pysam.AlignmentFile(cram_path, 'rc', reference_filename=args.reference_fasta)
        cram_filename = os.path.basename(cram_path)
        cram_filename_prefix = cram_filename.replace(".cram", "")
        total_reads = compute_total_reads(cram_path)

        record = {
            "SampleId": cram_filename_prefix,
            "TotalReads": total_reads,
        }
        smn1_nucleotide = nucleotide_count(cram, SMN1_C840_POSITION_1BASED - 1, SMN1_C840_POSITION_1BASED)
        record["SMN1_C_count"] = smn1_nucleotide["C"]
        record["SMN1_Total_count"] = sum(smn1_nucleotide.values())
        smn2_nucleotide = nucleotide_count(cram, SMN2_C840_POSITION_1BASED - 1, SMN2_C840_POSITION_1BASED)
        record["SMN2_C_count"] = smn2_nucleotide["C"]
        record["SMN2_Total_count"] = sum(smn2_nucleotide.values())
        record["SMN c.840: Reads with C"] = record["SMN1_C_count"]+record["SMN2_C_count"]
        record["SMN c.840: Total Reads"] = record["SMN1_Total_count"]+record["SMN2_Total_count"]
        if args.verbose:
            print(f"{total_reads:15,d} total reads in {cram_path}")
        output_records.append(record)
        cram.close()
    """
    scatter plot of the data generated in tsv file
    """
    df = pd.DataFrame(output_records)
    df.to_csv(args.output_tsv, sep='\t', index=False)
    print("No of cram files: ", len(output_records), " output in", args.output_tsv)
    df = pd.read_csv("c_count.tsv", sep='\t')
    SMN_C_count = df['SMN c.840: Reads with C']
    SMN_total = df['SMN c.840: Total Reads']
    x = list(SMN_C_count)
    y = list(SMN_total)
    plt.scatter(x=x,y=y,color = 'b')
    plt.grid(True)
    plt.xlabel('SMN c.840: Reads with C')
    plt.ylabel('SMN c.840: Total Reads')
    plt.savefig('C_vs_total.png')

if __name__ == "__main__":
    main()