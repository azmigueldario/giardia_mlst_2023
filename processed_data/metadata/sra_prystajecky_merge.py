#!/usr/bin/env python

# dependencies
import pandas as pd
import os
import argparse

# get script directory for default paths resolution
script_dir = os.path.dirname(os.path.abspath(__file__))

# main function
def merge_metadata_giardia(sra_metadata_path, paper_metadata_path, output_path):
    """
    Loads two metadata files, sanitizes column names, performs a left merge, 
    and saves the result to a new CSV file.
    """
    
    # 1. Load Data
    print(f"Loading Giardia metadata from: {sra_metadata_path}")
    try:
        sra_data = pd.read_csv(sra_metadata_path, sep="\t") 
    except FileNotFoundError:
        print(f"Error: Giardia metadata file not found at {sra_metadata_path}")
        return
        
    print(f"Loading paper metadata from: {paper_metadata_path}")
    try:
        paper_data = pd.read_csv(paper_metadata_path)
    except FileNotFoundError:
        print(f"Error: Paper metadata file not found at {paper_metadata_path}")
        return

    # Sanitize: strip whitespace, replace spaces/hyphens with underscores, and lowercase
    paper_data.columns = paper_data.columns.str.strip().str.replace(" ", "_").str.replace("-", "_")
    paper_data.columns = paper_data.columns.str.lower()
    
    # Keep first 6 columns (0:6) from paper data
    paper_data = paper_data.iloc[:, 0:6]
    
    # Left outer Merge
    print("Performing outer merge...")
    merged_data = pd.merge(
        sra_data, 
        paper_data, 
        left_on="sample_alias", right_on="ubc_id", 
        how="left",
        validate="many_to_one" )
    
    # drop redundant column
    merged_data = merged_data.drop(columns="ubc_id")
    
    # export Data to CSV
    print(f"Exporting merged data to: {output_path}")
    merged_data.to_csv(output_path, index=False)
    print("Merge complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge SRA metadata with Prystajecky paper metadata."
    )
    
    # Required arguments for input files
    parser.add_argument(
        "sra_metadata",
        nargs="?",
        default=os.path.join(script_dir, "Giardia_metadata_manual.tsv"),
        help="Absolute path to the Giardia metadata TSV file. Defaults to Giardia_metadata_manual.tsv."
    )
    parser.add_argument(
        "paper_metadata",
        nargs="?",
        default=os.path.join(script_dir, "prystajecky_2015_SuppData.csv"),
        help="Absolute path to the Prystajecky paper metadata CSV file. Defaults to prystajecky_2015_SuppData.csv"
    )
    
    # Argument for output file path with a default name
    parser.add_argument(
        "-o", "--output", 
        default=os.path.join(script_dir, "merged_SraPrystajecky.tsv"), 
        help="Name of the output TSV file. (Default: merged_SraPrystajecky.tsv)"
    )
    
    args = parser.parse_args()
    
    merge_metadata_giardia(
        args.sra_metadata, 
        args.paper_metadata, 
        args.output
    )
