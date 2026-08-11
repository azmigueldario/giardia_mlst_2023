#!/usr/bin/env python

####################################################################################
#                   DEPENDENCIES
####################################################################################
import os
import argparse
import sys
import logging
import pandas as pd
import numpy as np

# for testing
"""
sra_metadata_path = "/project/60006/mdprieto/giardia_mlst_2023/processed_data/metadata/sra_biosample_metadata_2025.csv"
paper_metadata_path = "/project/60006/mdprieto/giardia_mlst_2023/processed_data/metadata/prystajecky_2015_SuppData.csv"
output_path="/project/60006/mdprieto/giardia_mlst_2023/processed_data/metadata/merged_SraPrystajecky.csv"
colnames_drop=['secondary_', 'run_alias', '_title', 'md5', 'fastq_', 'datastore', 'create_date', 'version',
                'biosamplemodel', 'releasedate', 'lat_lon', 'isolation_source', 'projectid]
"""

####################################################################################
#                   Manual Fixes
####################################################################################

# define input and outputs
row_fixes = {
    'ATCC:50163': {'source': 'Cat - ATCC', 'location': 'PA,USA'},
    'ATCC:50803': {'source': 'Strain WB, clone 6 - ATCC', 'location': 'MD,USA'},
    'ATCC:30888': {'source': 'Human - ATCC', 'location': 'OR,USA'},
    'ATCC:50170': {'source': 'Animal - ATCC', 'location': 'WI,USA'},
    'Giardia AWB': {'source': 'Strain WB - ATCC', 'location': 'MD,USA'},
    'Giardia beaver': {'source': 'Beaver - donated'},
    'Giardia BGS': {'source': 'Strain GS - ATCC'}
}

####################################################################################
#               HELPER FUNCTIONS
####################################################################################

def parse_args():
    """
    Parse the arguments provided in the command-line
    """

    # get script directory for default paths resolution 
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # logging configuration
    logging.basicConfig(
        format="\n{levelname}:{message}",
        level=logging.DEBUG,
        style="{"
    )

    parser = argparse.ArgumentParser(description="Merge SRA metadata with Prystajecky paper metadata.")
    
    # Required arguments for input files
    parser.add_argument(
        "sra_metadata",
        nargs="?",
        default=None,
        help="Absolute path to the Giardia metadata TSV file. Defaults to Giardia_metadata_manual.tsv."
    )
    parser.add_argument(
        "paper_metadata",
        nargs="?",
        default=None,
        help="Absolute path to the Prystajecky paper metadata CSV file. Defaults to prystajecky_2015_SuppData.csv"
    )
    
    # Argument for output file path with a default name
    parser.add_argument(
        "-o", "--output", 
        default=None, 
        help="Name of the output TSV file. (Default: merged_SraPrystajecky.csv)"
    )

    default_drop_cols=['secondary_', 'run_alias', '_title', 'md5', 'fastq_', 'datastore', 'create_date', 'version',
                       'biosamplemodel', 'insert', 'sampletype', 'hash', 'lat_lon', 'projectid', 'submission']

    # list of columns to drop before exporting final dataset
    parser.add_argument(
        "-d", "--colnames_drop",
        nargs='*',
        type=str,
        default=default_drop_cols,
        help="List of columns to drop after merging (Default: {', '.join(default_drop_cols)})"
    )

    args = parser.parse_args()

    # assign sensible default paths, if not provided
    if args.sra_metadata is None:
        args.sra_metadata = os.path.join(script_dir, "sra_biosample_metadata_2025.csv")
    
    if args.paper_metadata is None:
        args.paper_metadata = os.path.join(script_dir, "prystajecky_2015_SuppData.csv")

    if args.output is None:
        args.output = os.path.join(script_dir, "merged_SraPrystajecky.csv")

    return args

def standardize_columns(df):
    """
    Sanitization for colnames: strip spaces/hyphens, replace slashes, make lowercase.
    """

    df.columns = (
        df.columns
            .str.replace(r'[\\/]', ' ', regex=True)
            .str.strip()
            .str.replace(r'[\s,-]+', '_', regex=True)
            .str.lower()
            .str.replace(r'_+', '_', regex=True)
            .str.strip('_')
        )
    
    return df

def apply_manual_fixes(df, row_fixes, col_name='samplename'):
    """
    Hard-coded strain fixes, requires input of a nested dictionary in the form:

    {'samplename': {'var1':'var1_value', 'var2':'var2_value'}

    The colnames must exist in the dataframe.
    """
    
    if col_name not in df.columns:
        logging.info(f"Columns required for apply_manual_fixes() not available in Dataframe.")
        return df
    
    for key, fixes in row_fixes.items():
        # filter and select rows matching colname key
        mask = df[col_name].str.contains(key, case=False, na=False)
        
        # loop and apply fixes if it matches 'mask'
        if 'source' in fixes:
            df.loc[mask, 'source'] = fixes['source']

        # pick if column was standardized
        if 'location' in fixes:
            target_col = 'province/state,_country' if 'province/state,_country' in df.columns else 'province_state_country'
            df.loc[mask, target_col] = fixes['location']
       
    return df

def resolve_column_information(df, source_cols, final_colname):
    """
    Parse a column by combining from different columns

    df: dataframe
    source_cols: list to define priority to complete the column
    final_colname: final column name to use for merged data
    """

    # define priority order list and presence in dataframe
    existing_source_cols = [column for column in source_cols if column in df.columns]

    if not existing_source_cols:
        logging.info(f"Columns required for resolve_column_information() in '{final_colname}' not available in Dataframe")
        return df
    
    # subfunction for cleaning column
    def clean_column(series):
        return (
            series
            .astype(str)
            .str.lower()
            .str.strip()
            .replace(['missing', 'not applicable', 'not collected', 'none', 'nan', ''], np.nan)
        )

    # initialize with 'isolation_source'
    consolidated = clean_column(df[existing_source_cols[0]])
    
    # update with following columns
    for col in existing_source_cols[1:]:
        consolidated = consolidated.combine_first(clean_column(df[col]))
    
    # define final column name
    df[final_colname] = (
        consolidated
        .fillna('NA')
        .str.replace('homo sapiens', 'human', case=False)
    )

    # drop all source columns but 'final_colname'
    cols_to_drop = [col for col in existing_source_cols if col != final_colname]  
    if cols_to_drop:
        df.drop(columns=cols_to_drop, inplace=True)

    return df

def parse_location(df):
    """
    Parse location independently of paper data and SRA data. 
    Then, merge location and clean unnecessary columns
    """

    # sanitize geo_location from SRA data
    if 'geo_loc_name' in df.columns:

        sra_location = (
            df['geo_loc_name']
                .str.replace(",", ":")
                .str.split(":", expand=True)
            )
        
        # save into individual columns if split string was available
        df['sra_country'] = sra_location[0].str.strip() if 0 in sra_location else None
        df['sra_city']    = sra_location[1].str.strip() if 1 in sra_location else None
        df['sra_region']  = sra_location[2].str.strip() if 2 in sra_location else None


    # sanitize Paper data
    if 'province_state_country' in df.columns:
        paper_location = (
            df['province_state_country']
            .str.rsplit(',', n=1, expand=True)
            )
        
        # save second column into region (or NA if empty)
        df['country_geoloc'] = paper_location[1].fillna(paper_location[0]).str.strip()
        df['region_geoloc']  = paper_location[0].where(paper_location[1].notna()).str.strip()
 
    #  new_var with finalized country information
    df['country_geoloc'] = df['country_geoloc'].combine_first(df['sra_country'])
    df['region_geoloc'] = df['region_geoloc'].combine_first(df['sra_region'] )
    df['city_geoloc'] = df['city_location'].combine_first(df['sra_city'])

    # clean obsolete columns
    df = df.drop(columns=['city_location', 'province_state_country', 'sra_country',
                          'sra_city', 'sra_region', 'geo_loc_name'], axis=1)
    
    return df

def filter_and_deduplicate(df, colnames_drop):
    """Final cleaning: drop columns, handle missing data, and resolve hybrid assemblies."""

    # Remove WB strains and redundant ID column
    df = df[~df['samplename'].str.contains('wb', case=False, na=False)].copy()
    if 'ubc_id' in df.columns:
        df = df.drop(columns=['ubc_id'])

    # Drop user-specified patterns
    regex_pattern = '|'.join(colnames_drop)
    cols_to_remove = df.columns[df.columns.str.contains(regex_pattern, case=False, na=False)]
    df = df.drop(columns=cols_to_remove)

    # Threshold drop: remove columns with >90% missing data
    df = df.dropna(axis=1, thresh=int(0.10 * len(df))).fillna('NA')

    # Deduplicate hybrid assemblies
    df['platform'] = df['platform'].str.strip().str.capitalize()
        # Only mark as hybrid if there is more than one UNIQUE platform name per BioSample
    df['unique_platforms'] = df.groupby('biosample')['platform'].transform('nunique')
    df.loc[df['unique_platforms'] > 1, 'platform'] = 'hybrid'
    df = (df.drop_duplicates(subset='biosample', keep='last')
            .drop(columns=['unique_platforms']))
            
    # remove duplicated column, if they exists
    df = df.drop(columns=['source_y', 'source', 'host'], errors='ignore')

    # clean NA strings
    garbage_regex = r'(?i)^(' + '|'.join(['missing', 'not applicable', 'not collected', 'none', 'nan']) + r')$'
    df = df.replace(to_replace=garbage_regex, value=np.nan, regex=True)
    df = df.replace('', np.nan)
    
    return df

####################################################################################
#               MAIN FUNCTION
####################################################################################

def merge_metadata_giardia(sra_metadata_path, paper_metadata_path, 
                           output_path, colnames_drop):
    """
    Loads two metadata files, sanitizes column names, performs a left join, 
    and saves the result to a new CSV file.
    """
    
    logging.info(f"Loading datasets")
    try:
        sra_data = pd.read_csv(sra_metadata_path, sep=",")
        paper_data = pd.read_csv(paper_metadata_path)
    except FileNotFoundError:
        logging.error("Giardia metadata file(s) not found")
        raise ValueError("Necessary input. data not available. Cannot proceed.\n")

    # Sanitize: strip whitespace, replace spaces/hyphens with underscores, and lowercase
    paper_data = standardize_columns(paper_data)
    sra_data = standardize_columns(sra_data)

    # Keep first 6 columns (0:6) from paper data
    paper_data = paper_data.iloc[:, 0:6]
    
    # merge datasets
    merged_data = pd.merge(
        sra_data, 
        paper_data, 
        left_on="samplename", right_on="ubc_id", 
        how="left",
        validate="many_to_one"
        )
    logging.info("Successfully performed outer merge.")

    # Additional parsing: manual fixes and parsing of individual columns
      
    merged_data = apply_manual_fixes(merged_data,  row_fixes)

    location_colnames = ['isolation_source', 'host', 'source', 'env_material', 'strain']
    merged_data = resolve_column_information(merged_data, location_colnames, 'isolation_source')

    collection_date_colnames = ['collection_date', 'date_of_collection']
    merged_data = resolve_column_information(merged_data, collection_date_colnames, 'collection_date')
    
    merged_data = parse_location(merged_data)
    merged_data = filter_and_deduplicate(merged_data, colnames_drop)

    # export to csv
    logging.info(f"Exporting merged data to: {output_path}\n")
    merged_data.to_csv(output_path, index=False, na_rep='NA')

####################################################################################
#               Command line execution
####################################################################################

if __name__ == "__main__":
    
    args = parse_args()

    sys.exit(
        merge_metadata_giardia(
        args.sra_metadata, 
        args.paper_metadata, 
        args.output,
        args.colnames_drop)
    )
