#!/usr/bin/env python
# coding: utf-8

# dependencies 
import os
import sys
import argparse
import hdbscan
import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
from sklearn.manifold import TSNE
from sklearn.cluster import HDBSCAN


def parse_arguments():
    """Takes a square matrix of distances, embeds them into a 2D plane using t-SNE and
    clusters the datapoints using HDBSCAN. All parameters probably need to be tuned
    outside the main script (e.g. though jupyter notebook interactive visualization).
    
    Finally, we subsample from the resulting clusters. 
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', 
                        type=argparse.FileType('r'),
                        help="Path to input file with sourmash distances")
    parser.add_argument('-o', '--outfile',
                        default=sys.stdout, type=argparse.FileType('w'),
                        help="Output of subsampled samples (FOFN) from clusters")
    parser.add_argument('--random_seed', type=int, default=112233,
                        help="Random seed for reproducibility")
    parser.add_argument('--min_cluster_size', type=int, default=15,
                        help="Minimum cluster size for HDBSCAN")
    parser.add_argument('--min_samples', type=int, default=15,
                        help="Minimum samples for HDBSCAN")
    parser.add_argument('--subsample_fraction', type=float, default=0.5,
                        help="Fraction of samples to subsample from each cluster")

    args = parser.parse_args()
    
    return args

def data_import_cleaning(filepath):
    """Imports the distance matrix and sanitizes the labels.
    """

    # load the data    
    sour_data = pd.read_csv(filepath, sep=',')

    # sanitize the labels
    sour_data.columns = sour_data.columns.str.replace('.*/', '', regex=True)
    sour_data.index = sour_data.columns

    return sour_data

def tune_tsne(sour_data, random_seed, iterations=5000):
    """Tunes the hyperparameters of t-SNE and saves the final embeddings.
    """

    # variables to capture loss in KL_divergenec (set to infinite initially) and best perplexity
    best_loss = np.inf
    best_perplexity = None
    
    # perplexity range depends on data size
    # define range of perplexity values and number of iterations for tuning
    data_length = sour_data.shape[0]
    logging.info("--- Starting t-SNE Automated Tuning ---\n")
    perplexity_range = np.arange(10, data_length, 10)
    max_iter =  iterations

    for p in perplexity_range:
        try:
            # Initialize t-SNE model for tuning, use random_seed for reproducibility
            model = TSNE(
                n_components=2, 
                metric='precomputed', 
                init="random",
                perplexity=p, 
                max_iter=max_iter,
                random_state=random_seed, 
                learning_rate='auto'
            )
            
            # Fit the model (we only care about the resulting loss)
            model.fit(sour_data) 
            
            # set current loss in the loop
            current_loss = model.kl_divergence_
            logging.info(f"Perplexity (p={p}): Final KL-Divergence = {current_loss:.4f}\n")
            
            # evaluate current KL_loss parameter against the best one so far
            # if it is the best, save the perplexity value
            if current_loss < best_loss:
                best_loss = current_loss
                best_perplexity = p
                
        except ValueError as e:
            logging.warning(f"Skipping perplexity={p}: {e}")
            continue

    if best_perplexity is None:
        logging.error("Could not find a valid perplexity in the range.")
        # Fall back to a standard value if search fails
        best_perplexity = 50 

    logging.info(f"Selected Optimal Perplexity (p): {best_perplexity} (Minimum KL-Divergence: {best_loss:.4f})\n")
    
    # execute a final model with the chosen perplexity
    logging.info(f"Running final t-SNE embedding with p={best_perplexity} and max_iter={max_iter}\n")
    
    final_model = TSNE(
        n_components=2, 
        random_state=random_seed, 
        metric='precomputed', 
        init="random",
        perplexity=best_perplexity,
        max_iter=max_iter,
        learning_rate='auto'
    )

    # save embeddings to DataFrame with original index

    tsne_embeddings = final_model.fit_transform(sour_data)
    vector_mat = pd.DataFrame(tsne_embeddings, columns=['t-SNE-1', 't-SNE-2'])
    vector_mat.index = sour_data.index

    return vector_mat


def hdbscan_automated(embeddings: pd.DataFrame, min_cluster_size, min_samples, subsample_fraction, random_seed):
    """
    Automates HDBSCAN clustering. The parameters are not agnostic, so you must input 
    a minimum cluster size and minimum samples to evaluate.
    A reasonable default could be min_cluster_size=15 and min_samples=15.

    Args:
        embeddings (pd.DataFrame): The low-dimensional t-SNE embeddings.
        hyperparameter min_cluster_size (int): Minimum cluster size for HDBSCAN.
        hyperparameter min_samples (int): Minimum samples for HDBSCAN.

    Returns:
        HDBSCAN: The final fitted HDBSCAN model object with optimal parameters.

    """
    logging.info("--- Starting HDBSCAN Automated Tuning ---\n")
    
    vector_mat = embeddings.copy()
    X = embeddings.values
    
    logging.info(f"Selected Optimal HDBSCAN Parameters: min_cluster_size={min_cluster_size} and min_samples_cluster=${min_samples}\n")
    
    final_hdb = HDBSCAN(
        algorithm="auto", 
        cluster_selection_method="eom",
        min_cluster_size=min_cluster_size, 
        min_samples=min_samples,
    ).fit(X)
    
    vector_mat['hdbscan'] = final_hdb.labels_
    vector_mat['cluster'] = vector_mat['hdbscan'].replace(-1, np.nan)

    # stratified sampling, sample the same fraction from each cluster
    ## droplevel() removes index added by groupby
    strat_sample = (
        vector_mat.groupby('cluster')[['t-SNE-1', 't-SNE-2', 'hdbscan', 'cluster']]
        .apply(lambda x: x.sample(frac=subsample_fraction, random_state=random_seed))
        .droplevel(0) 
        )

    unclustered_df = vector_mat[vector_mat['hdbscan'] == -1] 
    final_set = pd.concat([unclustered_df, strat_sample], axis=0)
    
    return final_set

def main(filepath, outfile, random_seed, min_cluster_size, min_samples, subsample_fraction):
    """
    Main function to run t-SNE and HDBSCAN automated clustering and subsampling.
    Clustering does not take into account the initial representation of the data, 
    instead it relies on sampling the same fraction from each cluster after HDBSCAN clustering.

    Default percentage of subsampling is 50% from each cluster.
    """

    # parse arguments
    args = parse_arguments()
    
    # set up logging, format as time:errorlevel:message, output to stdout
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    
    # import and clean data
    sour_data = data_import_cleaning(filepath)
    
    # tune t-SNE and get embeddings
    tsne_embeddings = tune_tsne(sour_data, random_seed)
    
    # automated HDBSCAN clustering
    final_samples = hdbscan_automated(tsne_embeddings, min_cluster_size, min_samples,
                                      subsample_fraction, random_seed)
    
    # output the final subsampled set
    final_samples.to_csv(outfile, sep='\t', index=True)

    logging.info(f"Subsampled data with {final_samples.shape[0]} samples saved to {outfile.name}\n")

if __name__ == '__main__':

    args = parse_arguments()
    
    main(args.infile.name, 
         args.outfile, 
         args.random_seed, 
         args.min_cluster_size, 
         args.min_samples,
         args.subsample_fraction)