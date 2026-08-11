#!/usr/bin/env python3

# =====================================================================================
#           Dependencies and argparse
# =====================================================================================

import argparse
import logging
import sys
import pandas as pd

def parser_arguments():

    parser = argparse.ArgumentParser(description="""
        Identify highest mapping reference genomes and add them to a master sheet.
        Mark those with low amount of combined mapped reads (<70%).
        
        Example:
        python3 map_giardia_assemblages.py --abundance file.csv --master master_metadata.csv --output path/updated_master.csv

        """)
    parser.add_argument("-a", "--abundance", nargs='+',
                        required=True, help="Path to the mapping summary CSV/TSV (the wide format file).")
    parser.add_argument("-m", "--master", required=True, help="Path to the master data sheet CSV.")
    parser.add_argument("-o", "--output", required=True, help="Path to save the updated master CSV.")
    parser.add_argument("--sample-col", default="Sample", help="Name of the column in the master sheet containing the SRR/ERR IDs (default: 'Sample').")
    
    args = parser.parse_args()

    return args

# =====================================================================================
#           Helper functions
# =====================================================================================


def import_process_abundance(abundance_files):
    """
    1. Read the abundance file (auto-detects if it's tab or comma separated).
    2. Process abundance data and create dictionary for mapping assigned assemblage.
    3. If 70% or more of the reads were not assigned, assign assemblage unknown
    """
    
    logging.info(f"Starting processing of {len(abundance_files)} abundance file(s).")

    # Loop through input files (e.g., Illumina + ONT)
    df_list = []
    for file_path in abundance_files:
        logging.info(f"Reading abundance file: {file_path}")
        try:
            df = pd.read_csv(file_path, sep=None, engine="python")
            if "Genome" not in df.columns:
                logging.error(f"'Genome' column not found in {file_path}.")
                sys.exit(1)
            
            # Set Genome as index so they align perfectly during concatenation
            df = df.set_index("Genome")
            df_list.append(df)
        except Exception as e:
            logging.error(f"Error reading {file_path}: {e}")
            sys.exit(1)

    # Combine all abundance files horizontally
    logging.info("Combining abundance files and transposing...")
    abundance_df = pd.concat(df_list, axis=1)

    # transpose based on Genome column
    abundance_df = abundance_df.T

    logging.info("Filtering relative abundances and extracting sample IDs...")
    abundance_df = abundance_df[abundance_df.index.str.contains("Relative Abundance", case=False, na=False)]

    # Extract IDs matching (SRR|ERR[0-9]+), convert values to numeric (errors to NaN)
    extracted_ids = abundance_df.index.str.extract(r'((?:SRR|ERR)\d+)')[0]
    abundance_df.index = extracted_ids
    abundance_df = abundance_df.apply(pd.to_numeric, errors='coerce')

    logging.info("Calculating highest mapping reference genomes...")
    genomes_only = abundance_df.drop(columns=['unmapped'], errors='ignore')
    genomes_only = genomes_only.fillna(0)
    best_genomes = genomes_only.idxmax(axis=1)
    if 'unmapped' in abundance_df.columns:
        best_genomes.loc[abundance_df['unmapped'] > 70] = 'unknown'

    # Convert to a mapping dictionary: {'ERR311174': 'GCF_000002435.2_WB_genomic', ...}
    mapping_dictionary = best_genomes.to_dict()

    return mapping_dictionary

def import_process_metadata(metadata_sheet, sample_indicator, output_path, mapping_dictionary):
    """
    Load master metadata csv, process it, and modify according to results of abundance.

    Save updated metadata file
    """

    logging.info(f"Reading master metadata sheet: {metadata_sheet}")
    try:
        master_df = pd.read_csv(metadata_sheet)
    except Exception as e:
        logging.error(f"Error reading {metadata_sheet}: {e}")
        sys.exit(1)

    if sample_indicator not in master_df.columns:
        logging.error(f"Sample column '{sample_indicator}' not found in the master sheet.")
        sys.exit(1)

    logging.info("Mapping highest reference genomes to metadata...")
    master_df['giardia_assemblage_mapping'] = master_df[sample_indicator].map(mapping_dictionary)

    logging.info("Assigning clean assemblage labels (A, B, or other)...")
    def map_assemblage(genome_str):
        if pd.isna(genome_str) or genome_str == 'unknown':
            return 'unknown'
        elif 'GCF_000002435.2' in genome_str:
            return 'assemblage A'
        elif 'GCA_000498735.1' in genome_str:
            return 'assemblage B'
        return 'other'
    
    master_df['assigned_assemblage'] = master_df['giardia_assemblage_mapping'].apply(map_assemblage)

    # Save output
    logging.info(f"Saving updated metadata to: {output_path}")
    master_df.to_csv(output_path, index=False)
    
    logging.info(f"Success! Processed {len(mapping_dictionary)} samples total.")
    
def main(abundance_data, metadata_sheet, sample_indicator, output_path):

    abundance_mapping = import_process_abundance(abundance_data)

    import_process_metadata(metadata_sheet, sample_indicator,
                            output_path, abundance_mapping)

    return 0

if __name__ == "__main__":
    args = parser_arguments()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    sys.exit(
        main(args.abundance,
             args.master, 
             args.sample_col,
             args.output)
    )
