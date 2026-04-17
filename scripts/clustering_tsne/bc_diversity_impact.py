#!/usr/bin/env python
# coding: utf-8

# ---------------------------------------------------------------------------------------------------
#       Preparations
# ---------------------------------------------------------------------------------------------------

import warnings
import umap
import re
import argparse
import logging
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

#  SILENCE WARNINGS 
warnings.filterwarnings("ignore", message="using precomputed metric")
warnings.filterwarnings("ignore", message="n_jobs value 1 overridden")
warnings.filterwarnings("ignore", message="The figure layout has changed to tight")

# SET LOGGER
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ARGUMENT PARSER
def parse_args():
    """Handle command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Downsample genomic subclusters using Min-Max diversity sampling and random sampling."
    )
    
    # Input/Output
    parser.add_argument("-d", "--distance_matrix", 
                        required=True, 
                        help="Path to Sourmash/Jaccard distance matrix (CSV)")
    parser.add_argument("-m", "--metadata", 
                        required=True, help="Path to metadata/vector_mat (CSV)")
    parser.add_argument("-o", "--output_prefix", 
                        default="giardia_subset", help="Output filename prefix (default: giardia_subset)")

    # Fasta Path Logic
    parser.add_argument("--fasta_dir", 
                        required=True, 
                        help="Base directory where FASTA files are stored")
    parser.add_argument("--fasta_ext", 
                        default=".fasta", 
                        help="Extension of fasta files (default: .fasta)")

    # Parameters
    parser.add_argument("--keep", 
                        type=int, default=20, 
                        help="Number of samples to keep from the target subcluster (default: 20)")
    parser.add_argument("--threshold_agglomerative",
                        default=0.05,
                        type=float, 
                        help="Distance. threshold for agglomerative clustering (Default: 0.05)")
    parser.add_argument("--threshold_umap",
                        required=True, 
                        type=float, 
                        help="UMAP-2 threshold for identifying the endemic lineage (defined interactively)")
    parser.add_argument("--neighbors_umap",
                        required=True, 
                        type=float, 
                        help="UMAP-2 influential neighbors (defined interactively)")
    parser.add_argument("--seed", type=int, 
                        default=1113, help="Random seed for reproducibility (default: 1113)")
    
    return parser.parse_args()

# ---------------------------------------------------------------------------------------------------
#       Helper functions
# ---------------------------------------------------------------------------------------------------

def run_agglomerative(dist_matrix, value_threshold):
    """Performs Hierarchical Clustering on a precomputed distance matrix."""
    clusterer = AgglomerativeClustering(
        n_clusters=None, 
        distance_threshold=value_threshold,
        metric='precomputed', 
        linkage='average'
    )
    return clusterer.fit(dist_matrix)

def run_umap_projection(dist_matrix, n_neighbors, seed):
    """Reduces high-dimensional distance data to 2D for visualization/labeling."""
    reducer = umap.UMAP(
        n_neighbors=n_neighbors, 
        min_dist=0.1, 
        metric='precomputed', 
        random_state=seed
    )
    return reducer.fit_transform(dist_matrix)

def build_metadata_lookup(dataframe, cluster_labels):
    """
    Creates a master lookup dictionary where keys are IDs (Run or BioSample) 
    and values are the full metadata row + cluster ID.
    """
    # Attach cluster IDs to the metadata dataframe first
    dataframe = dataframe.copy()
    dataframe['cluster_id'] = cluster_labels
    
    lookup = {}
    
    for _, row in dataframe.iterrows():
        # Create a clean metadata packet for this sample
        metadata_packet = {
            'isolate': str(row.get('isolate', 'Unknown')).replace("/", "-"),
            'host': str(row.get('isolation_source', 'Unknown')),
            'region': str(row.get('region_geoloc', 'Unknown')).lower(),
            'cluster': row['cluster_id'],
            'collection_date': str(row.get('collection_date', 'Unknown'))
        }
        
        # Map this packet to BOTH the Run and the BioSample
        for identifier in [str(row.get('run')), str(row.get('biosample'))]:
            if identifier and identifier != 'nan':
                lookup[identifier] = metadata_packet
                
    return lookup

def label_subclusters(row, threshold):
    """Dynamic labeling based on UMAP-2 threshold."""
    raw_region = str(row['Region']).lower()
    is_bc = bool(re.search(r'\bbc\b', raw_region))
    cluster_id = str(row.get('cluster_id', ''))
    
    if cluster_id == "2" and is_bc and row['UMAP-2'] > threshold:
        return 'BC_Only'
    return 'BC Samples' if is_bc else 'Global Context'

def sample_min_max(df, dist_matrix, label, n_samples):
    subset_ids = df[df['BC_source'] == label].index.tolist()
    if len(subset_ids) <= n_samples: return subset_ids

    # Start with a random sample
    selected = [subset_ids[0]]
    remaining = subset_ids[1:]
    while len(selected) < n_samples:
        current_distances = dist_matrix.loc[remaining, selected]
        next_point = current_distances.min(axis=1).idxmax()
        selected.append(next_point)
        remaining.remove(next_point)
    return selected


# ---------------------------------------------------------------------------------------------------
#       Helper functions
# ---------------------------------------------------------------------------------------------------


def clustering_and_subsampling(distance_matrix, metadata, agglo_threshold, umap_neighbors, 
                               seed_value, umap_threshold, keep, output_minmax, output_random):
    
    logger.info("fImporting primary datasets: {distance_matrix} and {metadata}")
    # Import and clean identifiers in distance matrix
    distance_df = pd.read_csv(distance_matrix, sep=',')
    distance_df.columns = distance_df.columns.str.replace('.*/', '', regex=True)
    distance_df.index = distance_df.columns

    metadata_df = pd.read_csv(metadata, sep=',')

    agglomerative_df = run_agglomerative(distance_df, agglo_threshold)
    logger.info("Building metadata-identifier lookup maps")
    identifier_map = build_metadata_lookup(metadata_df, agglomerative_df.labels_)

    # Create results dataframe and add some relevant metadata
    results_df = pd.DataFrame(index=distance_df.index)
    results_df['cluster_id']    = results_df.index.map(lambda x: identifier_map.get(x, {}).get('cluster_id', -1))
    results_df['Region']        = results_df.index.map(lambda x: identifier_map.get(x, {}).get('region', 'unknown'))
    results_df['Host']          = results_df.index.map(lambda x: identifier_map.get(x, {}).get('host', 'unknown'))

    # Verify completeness
    unmapped_count = (results_df['cluster_id'] == -1).sum()
    if unmapped_count > 0:
        logger.warning(f"{unmapped_count} samples in matrix could not be mapped to metadata!")

    # Add UMAP coordinates to results data
    embedding = run_umap_projection(distance_df, umap_neighbors, seed_value)
    results_df['UMAP-1'] = embedding[:, 0]
    results_df['UMAP-2'] = embedding[:, 1]

    # Identify subcluster with only BC samples
    results_df['BC_source'] = results_df.apply(lambda r: label_subclusters(r, umap_threshold), axis=1)

    # Min-max sampling and exporting
    logger.info(f"Performing Min-Max diversity sampling (target n: {keep})")
    other_ids = results_df[results_df['BC_source'] != 'BC_Only'].index.tolist()
    final_ids_minmax = sample_min_max(results_df, distance_df, 'BC_Only', keep) + other_ids
    logger.info(f"Assembling final dataset and exporting to {output_minmax}")
    final_output = results_df.loc[final_ids_minmax].copy()

    # Random sampling the subcluster
    logger.info(f"Performing random sampling (target n: {keep}) and saving to {output_random}")
    random_list = (results_df[results_df['BC_source'] == 'BC_Only']
                        .sample(n=20, random_state=seed_value)
                        .index
                        .to_list()
                        )
    final_ids_random = random_list + other_ids


def main():
    args = parse_args()
    clustering_and_subsampling(args.distance_matrix, args.metadata, args.threshold_agglomerative, args.neighbors_umap,
                               args.seed, args.umap_threshold)

if __name__ == "__main__":
    main()
