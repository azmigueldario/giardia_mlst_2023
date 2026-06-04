#!/usr/bin/env python
# coding: utf-8

# ---------------------------------------------------------------------------------------------------
#       Preparations
# ---------------------------------------------------------------------------------------------------

import os
import re
import umap
import random
import logging
import warnings
import argparse
import pandas as pd
from pathlib import Path
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
    parser.add_argument("--fasta_extension",
                        nargs='+', default=['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz'],
                        help="Extension of fasta files (default: .fasta)")

    # Parameters
    parser.add_argument("--agglo_cluster_value",
                        type=int, default=2,
                        help="Chosen cluster value to subsample")
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
    parser.add_argument("--umap_mode", choices=['above', 'below'], default='above',
                        help="Filter samples 'above' or 'below' the UMAP-2 threshold")
    parser.add_argument("--neighbors_umap",
                        required=True,
                        type=int,
                        help="UMAP-2 influential neighbors (defined interactively)")
    parser.add_argument("--seed", type=int,
                        default=1113, help="Random seed for reproducibility (default: 1113)")

    return parser.parse_args()

# ---------------------------------------------------------------------------------------------------
#       Helper functions
# ---------------------------------------------------------------------------------------------------

def clean_id(name):
    """Helper to strip paths and all genomic extensions including .gz"""
    if not isinstance(name, str):
        return str(name)
    # Remove path
    name = os.path.basename(name)
    # Remove extensions: matches .fasta, .fa, .fna, .fastq, .fq and optional .gz
    return re.sub(r'\.(fasta|fa|fna|fastq|fq)(\.gz)?$', '', name, flags=re.IGNORECASE)

def map_fasta_folder(fasta_dir, fasta_extensions):
    """Create lookup dictionary to identify paths for every accession"""

    all_paths = []
    base_path = Path(fasta_dir)

    # Check every extension in the base path. Map basename to absolute path.
    for ext in fasta_extensions:
        all_paths.extend(list(base_path.rglob(f"*{ext}")))

    logger.info(f"Indexed {len(all_paths)} unique files.")
    return all_paths

def build_metadata_lookup(input_dataframe, all_fasta_paths):
    """
    Matches samples to files by checking if ID is a SUBSTRING of the filename.
    Also adds important metadata labels
    """
    dataframe = input_dataframe.copy()
    dataframe.columns = [c.lower() for c in dataframe.columns]

    parts_map = {clean_id(path.name).lower(): str(path.absolute()) for path in all_fasta_paths}

    for filepath in all_fasta_paths:
        name_clean = clean_id(filepath.name).lower()

        for part in re.split(r'[_.-]', name_clean):
            if part and part not in parts_map:
                parts_map[part] = str(filepath.absolute())

    lookup = {}

    for _, row in dataframe.iterrows():
        run_id = str(row.get('run', row.get('run_id', ''))).lower()
        bio_id = str(row.get('biosample', '')).lower()


        metadata_packet = {
            'isolate': str(row.get('isolate', 'Unknown')).replace("/", "-"),
            'host': str(row.get('isolation_source', 'Unknown')),
            'region': str(row.get('region_geoloc', row.get('region', 'Unknown'))).lower(),
            'collection_date': str(row.get('collection_date', 'Unknown')),
            'fasta_path': parts_map.get(run_id, parts_map.get(bio_id, "NOT_FOUND"))
        }

        for identifier in [run_id, bio_id]:
            if identifier and identifier not in ['nan', '', 'none']:
                lookup[identifier] = metadata_packet
    return lookup

def label_subclusters(row, threshold, mode, cluster_id_value):
    """"
    Flexible filter to select samples in BC only subcluster
    """
    is_bc = bool(re.search(r'\bbc\b', str(row.get('Region', '')).lower()))
    cluster_match = str(row.get('cluster_id')) == str(cluster_id_value)

    # FLEXIBLE LOGIC: Check direction based on user input
    if mode == 'above':
        threshold_met = row['UMAP-2'] > threshold
    else:
        threshold_met = row['UMAP-2'] < threshold

    if cluster_match and is_bc and threshold_met:
        return 'BC_Only'
    return 'BC Samples' if is_bc else 'Global Context'

def sample_min_max(df, dist_matrix, label, n_samples, seed_value):
    subset_ids = df[df['BC_source'] == label].index.tolist()
    if len(subset_ids) <= n_samples:
        return subset_ids

    # Start with a random sample
    rng = random.Random(seed_value)
    first_pick = rng.choice(subset_ids)
    selected    = [first_pick]
    remaining   = [s for s in subset_ids if s != first_pick]
    while len(selected) < n_samples:
        current_distances = dist_matrix.loc[remaining, selected]
        next_point = current_distances.min(axis=1).idxmax()
        selected.append(next_point)
        remaining.remove(next_point)
    return selected

def map_compound_id(compound_id, identifier_map):
        compound_id = str(compound_id).lower()

        if compound_id in identifier_map:
            return identifier_map[compound_id]

        parts = re.split(r'[_.-]', compound_id)
        for part in parts:
            if part in identifier_map:
                return identifier_map[part]

        return {}

def generate_samplesheet(ids, lookup_map, output_name):
    """Produces the 'sample_id,contig' CSV, handling compound IDs."""
    rows = []
    for given_id in ids:

        # Try whole ID and then split for compound IDs
        entry = lookup_map.get(given_id.lower())
        if not entry:
            parts = re.split(r'[_.-]', str(given_id).lower())
            for part in parts:
                if part in lookup_map:
                    entry = lookup_map[part]
                    break

        # If still not found, default to an empty dict to avoid crashes
        entry = entry if entry else {}
        fasta = entry.get('fasta_path', 'NOT_FOUND')

        rows.append({'sample_id': given_id, 'contig': fasta})

    pd.DataFrame(rows).to_csv(output_name, index=False, header = ["sample", "contig"])
    logger.info(f"Samplesheet created: {output_name}")

# ---------------------------------------------------------------------------------------------------
#       Clustering helper functions
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
    safe_neighbors = max(2, min(n_neighbors, len(dist_matrix) - 1))

    reducer = umap.UMAP(
        n_neighbors=safe_neighbors,
        min_dist=0.1,
        metric='precomputed',
        random_state=seed
    )
    return reducer.fit_transform(dist_matrix)

# ---------------------------------------------------------------------------------------------------
#       Main functions
# ---------------------------------------------------------------------------------------------------

def clustering_and_subsampling(distance_matrix,
                               metadata,
                               agglo_threshold,
                               umap_neighbors,
                               umap_threshold,
                               umap_mode,
                               seed_value,
                               keep,
                               agglo_cluster_value,
                               fasta_dir,
                               fasta_extensions,
                               output_prefix):

    logger.info(f"Importing primary datasets: {distance_matrix} and {metadata}")
    # Import and clean identifiers in distance matrix
    distance_df = pd.read_csv(distance_matrix, sep=',')
    distance_df.columns = [clean_id(col) for col in distance_df.columns]
    distance_df.index = distance_df.columns

    metadata_df = pd.read_csv(metadata, sep=',')

    logger.info("Building metadata-identifier lookup maps")
    all_paths = map_fasta_folder(fasta_dir, fasta_extensions)
    agglomerative_df = run_agglomerative(distance_df, agglo_threshold)
    _cluster_mapping = dict(zip(distance_df.index, agglomerative_df.labels_))
    identifier_map = build_metadata_lookup(metadata_df, all_paths)

    # Create results dataframe and add some relevant metadata
    results_df = pd.DataFrame(index=distance_df.index)
    results_df['cluster_id'] = agglomerative_df.labels_

    results_df['metadata_packet'] = results_df.index.map(lambda x: map_compound_id(x, identifier_map))

    results_df['Region'] = results_df['metadata_packet'].apply(lambda x: x.get('region', 'unknown'))
    results_df['Host']   = results_df['metadata_packet'].apply(lambda x: x.get('host', 'unknown'))

    # Verify completeness
    unmapped_count = (results_df['cluster_id'] == -1).sum()
    if unmapped_count > 0:
        logger.warning(f"{unmapped_count} samples in matrix could not be mapped to metadata!")

    # Add UMAP coordinates to results data
    embedding = run_umap_projection(distance_df, umap_neighbors, seed_value)
    results_df['UMAP-1'] = embedding[:, 0]
    results_df['UMAP-2'] = embedding[:, 1]

    # Identify subcluster with only BC samples
    results_df['BC_source'] = results_df.apply(
        lambda row: label_subclusters(
            row,
            umap_threshold,
            mode=umap_mode,
            cluster_id_value=agglo_cluster_value
            ),
        axis=1
        )
    other_ids = results_df[results_df['BC_source'] != 'BC_Only'].index.tolist()
    bc_only_pool = results_df[results_df['BC_source'] == 'BC_Only']
    n_to_keep = min(keep, len(bc_only_pool))

    # debugger
    counts = results_df['BC_source'].value_counts()
    logger.info(f"Subsampling summary:\n{counts}")
    counts_cluster = results_df['cluster_id'].value_counts()
    logger.info(f"Subsampling summary:\n{counts_cluster}")

    bc_count = results_df['Region'].str.contains('bc').sum()
    logger.info(f"Successfully mapped {bc_count} samples to BC region.")


    # Check counts for each gate
    logger.info(f"Gate 1 (Region is BC): {results_df['Region'].str.contains('bc', case=False).sum()}")
    logger.info(f"Gate 2 (Cluster Match): {(results_df['cluster_id'].astype(str) == str(agglo_cluster_value)).sum()}")

    # Min-max sampling and exporting
    logger.info(f"Performing Min-Max diversity sampling (target n: {n_to_keep})")
    final_ids_minmax = sample_min_max(results_df, distance_df, 'BC_Only', n_to_keep, seed_value) + other_ids

    # Random sampling the subcluster
    logger.info(f"Performing random sampling (target n: {n_to_keep}) and saving to {output_prefix}")
    bc_random = bc_only_pool.sample(n=n_to_keep, random_state=seed_value).index.tolist()
    final_ids_random = bc_random + other_ids

    # Generate samplesheets
    generate_samplesheet(final_ids_minmax, identifier_map, f"{output_prefix}_minmax_samplesheet.csv")
    generate_samplesheet(final_ids_random, identifier_map, f"{output_prefix}_random_samplesheet.csv")

# ---------------------------------------------------------------------------------------------------
#       Run from the command line
# ---------------------------------------------------------------------------------------------------

def main():
    args = parse_args()
    clustering_and_subsampling(
            args.distance_matrix,
            args.metadata,
            args.threshold_agglomerative,
            args.neighbors_umap,
            args.threshold_umap,
            args.umap_mode,
            args.seed,
            args.keep,
            args.agglo_cluster_value,
            args.fasta_dir,
            args.fasta_extension,
            args.output_prefix
        )

if __name__ == "__main__":
    main()
