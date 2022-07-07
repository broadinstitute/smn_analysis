#!/usr/env python3

"""This script counts the number of reads that overlap exons of SMN1 and SMN2 in the input CRAM file."""

import argparse
import os
import pysam
import pandas as pd
import matplotlib.pyplot as plt


EXON_COORDINATES = [
    ("chr5:70925087-70925184","chr5:70049669-70049766"),
    ("chr5:70938839-70938910","chr5:70063415-70063486"),
    ("chr5:70941389-70941508","chr5:70065965-70066084"),
    ("chr5:70942358-70942558","chr5:70066934-70067134"),
    ("chr5:70942718-70942870","chr5:70067294-70067446"),
    ("chr5:70944658-70944753","chr5:70069235-70069330"),
    ("chr5:70946066-70946176","chr5:70070641-70070751"),
    ("chr5:70951941-70951994","chr5:70076521-70076574"),
    ("chr5:70952439-70953015","chr5:70077019-70077595"),
]

EXON_LABELS = {
    0: "1",
    1: "2a",
    2: "2b",
    3: "3",
    4: "4",
    5: "5",
    6: "6",
    7: "7",
    8: "8",
}


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
    os.system(f"samtools idxstat {cram_path} > {idxstat_filename}")
    idxstat_df = pd.read_table(idxstat_filename, names=[
        "chrom", "chrom_length", "mapped_read_counts", "unmapped_read_counts",
    ])
    os.remove(idxstat_filename)

    total_reads = sum(idxstat_df.mapped_read_counts) + sum(idxstat_df.unmapped_read_counts)

    return total_reads


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
            "exons 1+2a+2b+3+4+5+6": 0,
            "exon 7+8 ": 0

        }
        for exon_count, (SMN1_exon_region, SMN2_exon_region) in enumerate(EXON_COORDINATES):
            exon_label = EXON_LABELS[exon_count]
            SMN1_exon_read_count = cram.count(region= SMN1_exon_region)
            SMN2_exon_read_count = cram.count(region= SMN2_exon_region)
            total_exon = SMN1_exon_read_count + SMN2_exon_read_count
            record[f"Exon{exon_label}Reads"] = total_exon
            if exon_count < 7:
                record["exons 1+2a+2b+3+4+5+6"] += total_exon

            else:
                record["exon 7+8 "] += total_exon
        if args.verbose:
           print(f"{total_reads:15,d} total reads in {cram_path}")
        output_records.append(record)
        cram.close()
    df = pd.DataFrame(output_records)
    df.to_csv(args.output_tsv, sep='\t', index=False)
    print("No of cram files: ", len(output_records), " output in", args.output_tsv)


if __name__ == "__main__":
    main()