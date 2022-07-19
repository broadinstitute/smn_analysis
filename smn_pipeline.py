"""Hail Batch pipeline for running the smn_analysis.py script on 1 or more samples."""

import hashlib
import os
import pandas as pd
from step_pipeline import pipeline, Backend, Localize, Delocalize
import sys

# GRCh38 coordinates for SMN1 and SMN2 regions
SMN1_REGION = "chr5:70925087-70953015"
SMN2_REGION = "chr5:70049669-70077595"
SMN1_AND_SMN2_REGION = "chr5:70049669-70953015"
SMN1_AND_SMN2_REGION_WITH_PADDING = "chr5:70048669-70954015"

DOCKER_IMAGE = "weisburd/smn_analysis@sha256:d363a0b05a475624e47aabd9c667911eb01655e762a801ea1f855a636978857c"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

OUTPUT_FILENAME_PREFIX = "smn-exon-counts"

bp = pipeline("SMN analysis", backend=Backend.HAIL_BATCH_SERVICE)

# parse args
p = bp.get_config_arg_parser()
p.add_argument("-s", "--sample-id", action="append",
               help="Optional sample id to process")
p.add_argument("-n", "--num-samples-to-process", type=int,
               help="Process only the 1st n samples from the sample table (useful for testing)")
p.add_argument("-o", "--output-dir", default="gs://gnomad-bw2/SMN_analysis_with_rakshya",
               help="Base output directory where to copy output .tsv files")
p.add_argument("sample_table",
               help="Path of tab-delimited table containing samples to process along with their cram paths and other "
                    "metadata. It should contain at least the columns: sample_id, cram_path, crai_path, genome_version")
args = bp.parse_args()

df = pd.read_table(args.sample_table)

if "individual_id" in df.columns:
    df = df.rename({"individual_id": "sample_id"}, axis=1)

# validate args
missing_columns = {"sample_id", "cram_path", "crai_path", "genome_version"} - set(df.columns)
if missing_columns:
    p.error(f"{missing_columns} columns are missing from {args.sample_table}. The table has columns: {set(df.columns)}")


print(f"Read {len(df)} rows from {args.sample_table}")
df = df[df.genome_version == 38]
df = df[~df.cram_path.isna() & ~df.crai_path.isna()]
df = df.drop_duplicates(subset=['cram_path', 'crai_path'])
df = df.sort_values("sample_id")
print(f"Found {len(df)} rows with a cram_path and genome_version == 38")

if args.sample_id:
    df = df[df.sample_id.isin(args.sample_id)]
    print(f"Found {len(df)} out of {len(args.sample_id)} requested sample ids")
    if len(df) == 0:
        sys.exit(1)
    elif len(df) < len(args.sample_id):
        print(f"WARNING: Couldn't find sample ids: {set(args.sample_id) - set(df.sample_id)}")

if args.num_samples_to_process:
    df = df.iloc[:args.num_samples_to_process]

analysis_id = ", ".join(df.sample_id)  # compute a hash of the sample ids being processed
analysis_id = hashlib.md5(analysis_id.encode('UTF-8')).hexdigest().upper()
analysis_id = analysis_id[:10]  # shorten

if not args.force:
    bp.precache_file_paths(os.path.join(args.output_dir, OUTPUT_FILENAME_PREFIX) + "*")

steps = []
print(f"Processing {len(df)} samples")
for _, row in df.iterrows():
    s1 = bp.new_step(
        f"SMN analysis: {row.sample_id}",
        arg_suffix="step1",
        image=DOCKER_IMAGE,
        cpu=0.25,
        memory="standard",
        #storage="25Gi",
        output_dir=args.output_dir,
        delocalize_by=Delocalize.COPY,
    )
    s1.switch_gcloud_auth_to_user_account()
    reference_fasta_input, reference_fasta_fai_input = s1.inputs(
        REFERENCE_FASTA_PATH, REFERENCE_FASTA_FAI_PATH, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    s1.command("cd /io/")
    s1.command("set -ex")

    cram_input = s1.input(row.cram_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET)
    crai_input = s1.input(row.crai_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET)  # Localize.HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET
    s1.command(f"ls -lh {cram_input}")

    local_cram_path = f"{row.sample_id}.cram"
    s1.command(f"time samtools view -T {reference_fasta_input} -C {cram_input} {SMN1_AND_SMN2_REGION_WITH_PADDING} > {local_cram_path}")
    s1.command(f"time samtools index {local_cram_path}")

    output_tsv_name = f"{OUTPUT_FILENAME_PREFIX}.{row.sample_id}.tsv"
    s1.command(
        f"time python3 -u /smn_analysis.py "
        f"-R {reference_fasta_input} "
        f"--output-tsv {output_tsv_name} "
        f"--verbose "
        f"{local_cram_path}"
    )
    s1.command("ls")
    s1.output(output_tsv_name, delocalize_by=Delocalize.COPY)
    steps.append(s1)

# combine tables into a single table
if len(steps) > 1:
    s2 = bp.new_step(
        f"Combine {len(steps)} tables",
        image=DOCKER_IMAGE,
        cpu=0.25,
        memory="standard",
        #storage="70Gi",
        output_dir=args.output_dir,
        delocalize_by=Delocalize.COPY,
        arg_suffix="step2",
    )
    s2.command("set -ex")

    combined_output_tsv_filename = f"combined_results.{analysis_id}.{len(df)}_samples.tsv"
    for i, step in enumerate(steps):
        s2.depends_on(step)
        tsv_input = s2.use_previous_step_outputs_as_inputs(step, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        if i == 0:
            s2.command(f"head -n 1 {tsv_input} > {combined_output_tsv_filename}")
        s2.command(f"tail -n +2 {tsv_input} >> {combined_output_tsv_filename}")

    s2.output(combined_output_tsv_filename, delocalize_by=Delocalize.COPY)
elif len(steps) == 1:
    combined_output_tsv_filename = output_tsv_name

bp.run()

# download results and merge with the sample metadata table
os.system(f"gsutil -m cp {os.path.join(args.output_dir, combined_output_tsv_filename)} .")
result_df = pd.read_table(combined_output_tsv_filename)
df = df[set(df.columns) - {"cram_path", "crai_path"}]
df_with_metadata = pd.merge(result_df, df, how="left", left_on="SampleId", right_on="sample_id")
df_with_metadata.to_csv(combined_output_tsv_filename, sep="\t", header=True, index=False)
print(f"Wrote {len(df_with_metadata)} rows to {combined_output_tsv_filename}")