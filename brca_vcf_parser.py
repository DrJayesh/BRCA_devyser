#!/usr/bin/env python3
"""
brca_vcf_parser.py

This script parses all VCF files in a given folder, extracts relevant fields,
and writes the results into Excel files inside an "outputExcel" subfolder that
is created in the parent directory of the folder containing the VCFs.

For each VCF:
1. Skip the metadata lines to locate the "#CHROM" header.
2. Read the remaining tab-separated table (header + variant rows).
3. Split the INFO column into separate columns for AF, CLINVARPAT, RNA_ACC, and CONSEQUENCES.
4. Add a "VariantID" column in the format "chr-pos-ref-alt".
5. Save the final DataFrame to an Excel file in ../outputExcel/ relative to the
   provided folder.
6. Add a "VCF_COUNT" column showing in how many VCF files the variant occurs as "<count>/<total files>".
7. Save the final DataFrame to an Excel file in ./<input_folder>/outputExcel/.
8. Automatically run GeneBe annotation on these Excel files and store the
   results in a ``gnb_anno`` folder alongside ``outputExcel``.


Usage:
    python brca_vcf_parser.py
"""

import os
import sys
import pandas as pd
import annotate_genebe

def parse_info(info_str):
    """
    Given the INFO string (e.g. "AF=0.9701;RS=rs1799943;CLINVARPAT=Benign;..."),
    extract the values for AF, CLINVARPAT, RNA_ACC, and CONSEQUENCES.
    Return a dict with those keys (values as strings, or empty if missing).
    """
    # Split into key=value entries
    entries = info_str.split(';')
    info_dict = { "AF": "", "CLINVARPAT": "", "RNA_ACC": "", "CONSEQUENCES": "" }
    for entry in entries:
        if '=' not in entry:
            continue
        key, val = entry.split('=', 1)
        if key in info_dict:
            info_dict[key] = val
    return info_dict

def build_variant_counts(vcf_files, folder):
    """Return a mapping of VariantID to the number of VCF files in which that
    variant appears."""
    counts = {}
    for filename in vcf_files:
        vcf_path = os.path.join(folder, filename)
        header_line = None
        try:
            with open(vcf_path, "r") as f:
                for i, line in enumerate(f):
                    if line.startswith("#CHROM"):
                        header_line = i
                        break
        except OSError:
            continue
        if header_line is None:
            continue
        try:
            df = pd.read_csv(vcf_path, sep='\t', skiprows=header_line, dtype=str)
        except Exception:
            continue
        if "#CHROM" in df.columns:
            df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
        required = ["CHROM", "POS", "REF", "ALT"]
        if not all(col in df.columns for col in required):
            continue
        variant_ids = (
            df["CHROM"].astype(str) + "-" + df["POS"].astype(str) + "-" +
            df["REF"].astype(str) + "-" + df["ALT"].astype(str)
        )
        for vid in set(variant_ids):
            counts[vid] = counts.get(vid, 0) + 1
    return counts

def process_vcf_file(vcf_path, output_dir, variant_counts, total_files):
    """
    Read a VCF at ``vcf_path``, parse it, and write to an Excel file in
    ``output_dir``. ``variant_counts`` should be a mapping of ``VariantID`` to
    the number of VCF files in the batch that contain that variant. ``total_files``
    is the total number of VCF files being processed.
    """
    # Derive a base name for the output Excel (e.g., "sample.vcf" -> "sample.xlsx")
    base_name = os.path.splitext(os.path.basename(vcf_path))[0]
    output_path = os.path.join(output_dir, base_name + ".xlsx")

    # Determine how many header/meta lines to skip by locating the line that
    # starts with "#CHROM". VCF files may contain a variable number of meta
    # lines, so relying on a fixed count (e.g. 25) can cause parsing errors.
    header_line = None
    try:
        with open(vcf_path, "r") as f:
            for i, line in enumerate(f):
                if line.startswith("#CHROM"):
                    header_line = i
                    break
    except OSError as e:
        print(f"Error opening {vcf_path}: {e}", file=sys.stderr)
        return

    if header_line is None:
        print(f"Error: Could not locate header line in {vcf_path}", file=sys.stderr)
        return

    # Read the VCF using the discovered header line.
    try:
        df = pd.read_csv(vcf_path, sep='\t', skiprows=header_line, dtype=str)
    except Exception as e:
        print(f"Error reading {vcf_path}: {e}", file=sys.stderr)
        return

    # If the "#CHROM" column exists, rename it to "CHROM" for convenience.
    if "#CHROM" in df.columns:
        df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

    # Ensure the required columns exist
    required_cols = ["CHROM", "POS", "REF", "ALT", "INFO"]
    for col in required_cols:
        if col not in df.columns:
            print(f"File {vcf_path} is missing required column '{col}'. Skipping.", file=sys.stderr)
            return

    # Parse the INFO column into separate columns
    info_expanded = df["INFO"].apply(lambda x: pd.Series(parse_info(x)))
    # Append the new columns to the DataFrame
    df = pd.concat([df, info_expanded], axis=1)

    # Convert the AF value from string to a numeric type
    df['AF'] = pd.to_numeric(df['AF'], errors='coerce')

    # Create the "VariantID" column as "chr-pos-ref-alt"
    df["VariantID"] = (
        df["CHROM"].astype(str) + "-" + df["POS"].astype(str) + "-" +
        df["REF"].astype(str) + "-" + df["ALT"].astype(str)
    )

    df["VCF_COUNT"] = df["VariantID"].apply(
        lambda vid: f"{variant_counts.get(vid, 0)}/{total_files}"
    )

    # Write to Excel
    try:
        with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
            df.to_excel(writer, index=False)
        print(f"Wrote parsed data to {output_path}")
    except Exception as e:
        print(f"Error writing Excel for {vcf_path}: {e}", file=sys.stderr)

def main():
    # Prompt the user for the folder containing VCF files
    folder = input("Enter the path to the folder containing VCF files: ").strip()
    if not os.path.isdir(folder):
        print(f"Error: '{folder}' is not a valid directory.", file=sys.stderr)
        sys.exit(1)

    # Create the output directory in the parent directory of the given folder
    parent_dir = os.path.dirname(os.path.abspath(folder))
    output_dir = os.path.join(parent_dir, "outputExcel")
    os.makedirs(output_dir, exist_ok=True)

    # Process each .vcf (or .VCF) file in the folder
    vcf_files = [f for f in os.listdir(folder) if f.lower().endswith(".vcf")]
    if not vcf_files:
        print("No VCF files found in the specified folder.")
        sys.exit(0)

    # First collect how many files contain each VariantID
    variant_counts = build_variant_counts(vcf_files, folder)

    for filename in vcf_files:
        vcf_path = os.path.join(folder, filename)
        process_vcf_file(vcf_path, output_dir, variant_counts, len(vcf_files))

    # Automatically annotate the produced Excel files using GeneBe
    annotate_genebe.annotate_folder(output_dir)

if __name__ == "__main__":
    main()
