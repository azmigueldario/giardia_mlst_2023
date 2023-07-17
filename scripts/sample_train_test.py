""" Script to produce n-fold cross validation datasets """

# requires numpy, pandas, os, argparse and sklearn

import os
import argparse
import pandas as pd
from sklearn.model_selection import train_test_split

# ---------------------------------- define parsing of arguments from command line
#def parse_arguments():
"""Reads arguments from the command line for my script."""

parser = argparse.ArgumentParser(description="Produces 'n' pairs of .csv files with testing and training datasets.")
parser.add_argument("-n", "--number_iterations", metavar='0,1,2,10,100...', 
                    default=10, required=False, type=int,
                    help='number of iterations for cross validation')
parser.add_argument("-o", "--output_dir", metavar='', required=True, type=str,
                    help='path [absolute/relative] to directory where output file(s) will be saved')
parser.add_argument("--test_proportion", default=0.2, type=float,
                    help="percentage of available data to be used for testing.")
parser.add_argument("--long_format", action="store_true",
                    help="outputs results as a single LONG format .csv file")
parser.add_argument("--separate", action="store_true",
                    help="outputs results as several .csv files. Produces wide format by default")
parser.add_argument("input_file", 
                    help="path [absolute/relative] to an input dataset in .csv format")

args = parser.parse_args()

# user defined inputs
directory = args.output_dir
n = args.number_iterations 

# sanity check for arguments
if args.test_proportion > 1 or args.test_proportion < 0:
    raise ValueError('--test_proportion must be between 0 and 1.')

# ---------------------------------  define function

def n_fold_split():
    """ This function takes a '.csv' data file as input. It divides it into 'training' and 'testing' datasets
        inside a for loop 'n' times. Every time a dataset is produced, the script saves it in the specified
        output directory.             
    """

    # read dataset
    df = pd.read_csv(f'{args.input_file}',
                     header=0)

    # create output_dir if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    

    if args.long_format:
        # to produce tests in long_format
        merged_long = pd.DataFrame()     # create empty dataframe

    
        for i in range(1, (n+1)):
        
            train, test = train_test_split (df, test_size=args.test_proportion)
            
                # assign a new column with selected values
            train = train.assign(set = f'set{i}', value = "train")
            test = test.assign(set = f'set{i}', value = "test")
            
                # concatenate the results of the iteration
            merged_long = pd.concat([merged_long, train, test], ignore_index=True)

            # export file
        merged_long.to_csv(f'{args.output_dir}/train_set_long.csv',
                           index=False)

    elif args.separate:
        # produces results in separate files
        for i in range(1, (n+1)):

            train, test = train_test_split (df,
                                            test_size=args.test_proportion)

            train.to_csv(f'{args.output_dir}/train_set_{i}.csv', index=False)
            test.to_csv(f'{args.output_dir}/test_set_{i}.csv', index=False)   

    else:
        # produces results in wide format by default
        merged_wide = df.copy()     # creates copy to avoid modifying primary one

        for i in range(1, 11):
            train, test = train_test_split (df, test_size=0.2) # split test/train
            train[f'set{i}'] = "train" # assign column names
            test[f'set{i}'] = "test"
                # merges vased on sample/contig columns, validates match 1:1 in rows or drops a column
            merged_wide= merged_wide.merge(pd.concat([train, test]), # concatenates train/test subsets
                                        on=["sample", "contig"],
                                        validate="1:1")

            # export file
        merged_wide.to_csv(f'{args.output_dir}/train_set_wide.csv',
                           index=False)


# run when script is called from CLI using 'python3 script.py'
if __name__ == '__main__':
    n_fold_split()
