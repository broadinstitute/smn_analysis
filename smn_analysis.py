#!/usr/env python3

"""This script counts the number of reads that overlap exons of SMN1 and SMN2 in the input CRAM file."""

import argparse
import os
import pysam
import pandas as pd
import matplotlib.pyplot as plt

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
    SMN_Nucleotide = {"T":0,"C":0,"G":0,"A":0,"N":0}
    for pileupcolumn in cram.pileup(contig="chr5", start=start, end=end, truncate=True):
        for pileupread in pileupcolumn.pileups:
            val = (pileupread.alignment.query_sequence[pileupread.query_position])
            if val.upper() in SMN_Nucleotide:
                SMN_Nucleotide[val] += 1
            else:
                raise Exception('Error')
    return SMN_Nucleotide

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
        SMN1_Nucleotide = nucleotide_count(cram,70951945,70951946)
        record["SMN1_C_count"] = SMN1_Nucleotide["C"]
        record["SMN1_Total_count"] = sum(SMN1_Nucleotide.values())
        SMN2_Nucleotide = nucleotide_count(cram,70076525,70076526)
        record["SMN2_C_count"] = SMN2_Nucleotide["C"]
        record["SMN2_Total_count"] = sum(SMN2_Nucleotide.values())
        record["SMN c.840:C"] = record["SMN1_C_count"]+record["SMN2_C_count"]
        record["SMN c.840:total"] = record["SMN1_Total_count"]+record["SMN2_Total_count"]
        if args.verbose:
            print(f"{total_reads:15,d} total reads in {cram_path}")
        output_records.append(record)
        cram.close()
    df = pd.DataFrame(output_records)
    df.to_csv(args.output_tsv, sep='\t', index=False)
    print("No of cram files: ", len(output_records), " output in", args.output_tsv)


if __name__ == "__main__":
    main()