#!/usr/env python3

"""Counts and prints the number of reads that overlap exon 7 of SMN1 and SMN2 in the input CRAM or BAM file. 
"""
import argparse
import os
import pysam
import pandas as pd


# parse arguments
EXON_COORDINATES = [("chr5:70925087-70925184","chr5:69345439-69345593"),
                    ("chr5:70938839-70938910","chr5:69359242-69359313"),
                    ("chr5:70941389-70941508","chr5:69361792-69361911"),
                    ("chr5:70942358-70942558","chr5:69362761-69362961"),
                    ("chr5:70942718-70942870","chr5:69363121-69363273"),
                    ("chr5:70944658-70944753","chr5:69365062-69365157"),
                    ("chr5:70946066-70946176","chr5:69366468-69366578"),
                    ("chr5:70951941-70951994","chr5:69372348-69372401"),
                    ("chr5:70952439-70953015","chr5:69372846-69373419")]

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
    parser.add_argument("-o", "--output-tsv", default="add-all-exon-count.tsv", help="Output tsv file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Whether to print extra stats and output")
    parser.add_argument("cram_paths", nargs="+", help="One or more CRAM file paths")
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
        #listing exon coordinates exon 1-8
        cram_filename = os.path.basename(cram_path)
        cram_filename_prefix = cram_filename.replace(".cram", "")
        #empty dictionary
        record = {}
        exon_count = 0
        #updating the records dictionary
        record.update({"SAMPLEID": cram_filename_prefix})
        #looping to read exon coordinates from the list
        for SMN1_exon_region, SMN2_exon_region in EXON_COORDINATES:
            SMN1_exon_read_count = cram.count(region= SMN1_exon_region)
            SMN2_exon_read_count = cram.count(region= SMN2_exon_region)
            exon_label = EXON_LABELS[exon_count]
            exon_count += 1
            record.update({
                f"SMN1 + SMN2 exon {exon_label} reads": SMN1_exon_read_count + SMN2_exon_read_count
            })

        total_reads = compute_total_reads(cram_path)
        record.update({
            "Total reads": total_reads,
        })

        if args.verbose:
            print(f"{total_reads:15,d} total reads in {cram_path}")

        output_records.append(record)
        cram.close()

    df = pd.DataFrame(output_records)
    df.to_csv(args.output_tsv, sep='\t', index=False)
    print("No of cram files: ", len(output_records), " output in", args.output_tsv)


if __name__ == "__main__":
    main()