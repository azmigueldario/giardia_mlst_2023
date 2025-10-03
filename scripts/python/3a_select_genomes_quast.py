#!/usr/bin/env python3

# Import dependencies 
import os
import sys
import argparse
import pandas as pd

# Define constants
GENOME_SIZE = 12000000

def parse_arguments():
    """Filters the results of assembly quality control in a .csv file and outputs
    a list of genomes in a text file (one file per line).
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', 
                        help="Path to input file", type=argparse.FileType('r'),   )
    parser.add_argument('-o', '--outfile', 
                        help="Output file",
                        default=sys.stdout, type=argparse.FileType('w'))
    parser.add_argument('--n50', type=int, default=30000,
                        help="Threshold value for N50")
    parser.add_argument('--contigs_n', type=int, default=1300,
                        help="Maximum value of contigs in assembly")
    parser.add_argument('--genome_size', type=int, default=GENOME_SIZE,
                        help="Length in basepairs of the reference genome.")

    args = parser.parse_args()
    
    return args

def parse_csv(filepath, output_file, n50, contigs_n, genome_size):
    """"Evaluates the input csv and selects the rows that meet the specified criteria.
    It returns a txt file with a list of filenames/paths.
    """
    data = pd.read_csv(filepath, sep='\t')
    
    condition_1 = data['# contigs'] < contigs_n
    condition_2 = data['N50'] > n50
    condition_3 = (data['Total length'] > genome_size*0.8) & (data['Total length'] < genome_size*1.2)
    
    filtered_df = data [condition_1 & condition_2 & condition_3]

        # output the ID column
    if output_file is not sys.stdout:
        (filtered_df
         .iloc[:, [0]]
         .to_csv(output_file, sep='\t', header=False, index=False)
         )
    else:
        print(filtered_df.iloc[:, [0]])
        print("Total samples:" + str(len(filtered_df.index)))
    
    
if __name__ == '__main__':
    args = parse_arguments()
    parse_csv(args.infile, args.outfile, args.n50, args.contigs_n, args.genome_size)