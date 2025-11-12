#!/usr/bin/env python

# dependencies
import pandas as pd
import os
import argparse

# get script directory for default paths resolution
script_dir = os.path.dirname(os.path.abspath(__file__))


# main function
def merge_metadata_giardia(sra_metadata_path, paper_metadata_path, 
                           output_path, colnames_drop):
    """
    Loads two metadata files, sanitizes column names, performs a left merge, 
    and saves the result to a new CSV file.
    """
    
    # 1. Load Data
    print(f"\nLoading Giardia metadata from: {sra_metadata_path}")
    try:
        sra_data = pd.read_csv(sra_metadata_path, sep="\t") 
    except FileNotFoundError:
        print(f"\nError: Giardia metadata file not found at {sra_metadata_path}")
        return
        
    print(f"\nLoading paper metadata from: {paper_metadata_path}")
    try:
        paper_data = pd.read_csv(paper_metadata_path)
    except FileNotFoundError:
        print(f"\nError: Paper metadata file not found at {paper_metadata_path}")
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
        validate="many_to_one"
        )
    
    # manual from NCBI BioSample descriptions
    merged_data.loc[merged_data['sample_alias'].str.contains('ATCC:50163', regex=True, case=False), ['source', 'province/state,_country']] = ['Cat - ATCC', 'PA,USA']
    merged_data.loc[merged_data['sample_alias'].str.contains('ATCC:50803', regex=True, case=False), ['source', 'province/state,_country']] = ['Strain WB, clone 6 - ATCC', 'MD,USA']
    merged_data.loc[merged_data['sample_alias'].str.contains('ATCC:30888', regex=True, case=False), ['source', 'province/state,_country']] = ['Human - ATCC', 'OR,USA']
    merged_data.loc[merged_data['sample_alias'].str.contains('ATCC:50170', regex=True, case=False), ['source', 'province/state,_country']] = ['Animal - ATCC', 'WI,USA']
    merged_data.loc[merged_data['sample_alias'] == 'Giardia beaver', ['source']] = ['Beaver - donated']
    merged_data.loc[merged_data['sample_alias'] == 'Giardia BGS', ['source']] = ['Strain GS - ATCC']
    merged_data.loc[merged_data['sample_alias'] == 'Giardia AWB', ['source', 'province/state,_country']] = ['Strain WB - ATCC', 'MD,USA']
      
    # drop redundant ID column
    merged_data = merged_data.drop(columns="ubc_id")

    # drop other columns unlikely to be used in downstream analysis
    print(f"\n\nRemoving unnecessary columns from merged data: {colnames_drop}")
    regex_drop_pattern = '|'.join(colnames_drop)
    columns_to_drop = merged_data.columns[merged_data.columns.str.contains(regex_drop_pattern, case=False, na=False)].tolist()
    merged_data = merged_data.drop(columns=columns_to_drop, axis=1)

    # explicit NA
    merged_data = merged_data.fillna('NA')
    
    # export Data to CSV
    print(f"\nExporting merged data to: {output_path}")
    merged_data.to_csv(output_path, index=False)
    print(f"\nMerge complete!")


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
        default=os.path.join(script_dir, "merged_SraPrystajecky.csv"), 
        help="Name of the output TSV file. (Default: merged_SraPrystajecky.csv)"
    )

    parser.add_argument(
        "-d", "--colnames_drop",
        nargs='*',
        type=str,
        default=['secondary_', 'run_alias', '_title', 'md5', 'fastq_'],
        help="List of columns to drop after merging (Default: secondary_  run_alias _title md5 bytes galaxy aspera fastq_ftp])"
    )
    
    args = parser.parse_args()
    
    merge_metadata_giardia(
        args.sra_metadata, 
        args.paper_metadata, 
        args.output,
        args.colnames_drop
    )
