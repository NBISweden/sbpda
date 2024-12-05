#!/usr/bin/env python
import argparse
import os
import random
import warnings
from datetime import datetime

import anndata as ad
import gdown
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import squidpy as sq
import torch 
from matplotlib import gridspec
from sklearn.preprocessing import MinMaxScaler

from nichecompass.models import NicheCompass
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                create_new_color_dict,
                                compute_communication_gp_network,
                                visualize_communication_gp_network,
                                extract_gp_dict_from_mebocost_ms_interactions,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                filter_and_combine_gp_dict_gps_v2,
                                generate_enriched_gp_info_plots)


 

# AnnData keys
counts_key = "counts"
adj_key = "spatial_connectivities"
gp_names_key = "nichecompass_gp_names"
active_gp_names_key = "nichecompass_active_gp_names"
gp_targets_mask_key = "nichecompass_gp_targets"
gp_targets_categories_mask_key = "nichecompass_gp_targets_categories"
gp_sources_mask_key = "nichecompass_gp_sources"
gp_sources_categories_mask_key = "nichecompass_gp_sources_categories"
latent_key = "nichecompass_latent"

# Architecture
cat_covariates_embeds_injection = ["gene_expr_decoder"]
cat_covariates_no_edges = [True]

# Analysis 
differential_gp_test_results_key = 'nichecompass_differential_gp_test_results'


def load_samplesheet(samplesheet):
    if not samplesheet or not os.path.exists(samplesheet):
        print("Error: samplesheet is missing or does not exist. Please provide a valid CSV file with path information.")
        sys.exit(1)
    else:
        print(f"Samplesheet found at: {samplesheet}")
        return samplesheet

def prepare_folders(base_path, data_path, tempDir, species):
    """
    Prepares the necessary folder structure for the project.
    This function creates the base, data, and temporary directories if they do not exist.
    It also sets up specific subdirectories for gene annotations, gene programs, artifacts,
    models, and figures. The function ensures that all paths are absolute and checks for
    the existence of the samplesheet.
    Parameters:
    base_path (str): The base directory path.
    data_path (str): The data directory path.
    tempDir (str): The temporary directory path.
    species (str): The species identifier used for naming specific files.
    samplesheet (str): The path to the samplesheet CSV file.
    Returns:
    tuple: A tuple containing the paths to the figure folder and model folder.
    Raises:
    SystemExit: If the samplesheet is missing or does not exist, or if there is an error
                creating the base directory.
    """

    
    # Use absolute paths
    base_path = os.path.abspath(base_path)
    data_path = os.path.abspath(data_path)
    tempDir = os.path.abspath(tempDir)

    print(f"Attempting to create or access base directory: {base_path}")

    try:
        os.makedirs(base_path, exist_ok=True)
        print(f"Base directory created or already exists: {base_path}")
    except Exception as e:
        print(f"Error creating base directory: {e}")
        sys.exit(1)

    # Check and create base, data, and temporary directories if they do not exist
    for path, label in zip([data_path, tempDir], ['Data directory', 'Temporary directory']):
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"{label} created at: {path}")
        else:
            print(f"{label} exists at: {path}")

    
    # Create folders for gene programs and gene annotations
    ga_data_folder_path = os.path.join(data_path, 'gene_annotations')
    gp_data_folder_path = os.path.join(data_path, 'gene_programs')
    artifacts_folder_path = os.path.join(tempDir)

    # Create all necessary folders if they do not exist
    os.makedirs(ga_data_folder_path, exist_ok=True)
    os.makedirs(gp_data_folder_path, exist_ok=True)
    os.makedirs(artifacts_folder_path, exist_ok=True)
    
    # Define paths for specific resources
    omnipath_lr_network_file_path = os.path.join(gp_data_folder_path, 'omnipath_lr_network.csv')
    nichenet_lr_network_file_path = os.path.join(gp_data_folder_path, f'nichenet_lr_network_v2_{species}.csv')
    nichenet_ligand_target_matrix_file_path = os.path.join(gp_data_folder_path, f'nichenet_ligand_target_matrix_v2_{species}.csv')
    mebocost_enzyme_sensor_interactions_folder_path = os.path.join(gp_data_folder_path, 'metabolite_enzyme_sensor_gps')
    # Check whether you need 'human_mouse_gene_orthologs.csv' file ################????????????????????
    gene_orthologs_mapping_file_path = os.path.join(ga_data_folder_path, 'human_mouse_gene_orthologs.csv')

    # Create model and figures directories within artifacts
    figure_folder_path = os.path.join(artifacts_folder_path, 'figures')
    model_folder_path = os.path.join(artifacts_folder_path, 'model')
    os.makedirs(figure_folder_path, exist_ok=True)
    os.makedirs(model_folder_path, exist_ok=True)
    
    print(f"Figure folder path: {figure_folder_path}")
    print(f"Model folder path: {model_folder_path}")

    return figure_folder_path, model_folder_path

def load_anndata(samplesheet_path):
    valid_files_cntr = 0
    missing_files_cntr = 0
    missing_files = []
    valid_files = []
    
    samples = pd.read_csv(samplesheet_path)
    for fl in range(samples.shape[0]):
        file_path = os.path.join(data_path, samples['adata_path'][fl])
        if os.path.exists(file_path):
            valid_files_cntr = valid_files_cntr + 1
            valid_files.append(file_path)
        else:
            missing_files_cntr = missing_files_cntr + 1
            missing_files.append(samples['adata_path'][fl])
    if missing_files_cntr > 0:
        print(f'There are {valid_files_cntr} / {samples.shape[0]} files exist. Please check the following samples and rerun.\n {missing_files}')
    else:
        print(f'All {valid_files_cntr} / {samples.shape[0]} files exist')

    for batch in valid_files:
        print(f"Loading batch {batch}...")

        adata_batch = sc.read_h5ad(batch)
        print("Computing spatial neighborhood graph...\n")
        # Compute (separate) spatial neighborhood graphs
        sq.gr.spatial_neighbors(adata_batch,
                                coord_type="generic",
                                spatial_key=spatial_key,
                                n_neighs=n_neighbors)
        
        # Make adjacency matrix symmetric
        adata_batch.obsp[adj_key] = (
            adata_batch.obsp[adj_key].maximum(
                adata_batch.obsp[adj_key].T))
        adata_batch_list.append(adata_batch)
    adata = ad.concat(adata_batch_list, join="inner")
    sc.write(os.path.join(model_folder_path, "adata_merged.h5ad"), adata)
    return adata

# Separate the data preparation and analysis steps:
    # [x] samplesheet
    # [x] folder preparation
    # [] loading the anndata objects
    # [] Combine spatial neighborhood graphs as disconnected components
    # [] Extract gene programs from gene-gene interaction networks
    # [] Filter and combine gene programs
    # [] Train the NicheCompass model
    # [] Analyze the trained model
    # [] Visualize the results
    # [] Save the results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare necessary folders for NicheCompass analysis.")
    parser.add_argument('--base_path', required=True, help="Base directory for the project")
    parser.add_argument('--data_path', required=True, help="Directory for storing gene annotations and programs data")
    parser.add_argument('--tempDir', required=True, help="Temporary directory for storing intermediate files")
    parser.add_argument('--species', default="human", help="Species for NicheCompass analysis (default: human)")
    parser.add_argument('--samplesheet', help="Sample sheet file")

    # Additional parameters
    parser.add_argument('--spatial_key', help="Spatial key for analysis")
    parser.add_argument('--n_neighbors', type=int, help="Number of neighbors")
    
    # AnnData keys
    parser.add_argument('--cat_covariates_keys', help="Categorical covariates keys")
    
    # Architecture
    parser.add_argument('--cat_covariates_embeds_nums', help="Categorical covariates embed numbers")
    parser.add_argument('--conv_layer_encoder', help="Convolution layer encoder")
    parser.add_argument('--active_gp_thresh_ratio', type=float, help="Active GP threshold ratio")
    
    # Trainer
    parser.add_argument('--n_epochs', type=int, help="Number of epochs")
    parser.add_argument('--n_epochs_all_gps', type=int, help="Number of epochs for all GPs")
    parser.add_argument('--lr', type=float, help="Learning rate")
    parser.add_argument('--lambda_edge_recon', type=float, help="Lambda edge reconstruction")
    parser.add_argument('--lambda_gene_expr_recon', type=float, help="Lambda gene expression reconstruction")
    parser.add_argument('--lambda_l1_masked', type=float, help="Lambda L1 masked")
    parser.add_argument('--lambda_l1_addon', type=float, help="Lambda L1 addon")
    parser.add_argument('--edge_batch_size', type=int, help="Edge batch size")
    parser.add_argument('--n_sampled_neighbors', type=int, help="Number of sampled neighbors")
    parser.add_argument('--use_cuda_if_available', type=bool, help="Use CUDA if available")
    
    # Analysis
    parser.add_argument('--cell_type_key', help="Cell type key")
    parser.add_argument('--latent_leiden_resolution', type=float, help="Latent Leiden resolution")
    parser.add_argument('--latent_cluster_key', help="Latent cluster key")
    parser.add_argument('--sample_key', help="Sample key")
    parser.add_argument('--spot_size', type=int, help="Spot size")
    
    args = parser.parse_args()

    # Run folder preparation and retrieve paths
    figure_folder_path, model_folder_path = prepare_folders(
        base_path=args.base_path,
        data_path=args.data_path,
        tempDir=args.tempDir,
        species=args.species
    )

    samplesheet_path = load_samplesheet(args.samplesheet)
    adata = load_anndata(samplesheet_path)
    print(f"adata is saved in {model_folder_path}")



    # Write the paths to files
    with open('figure_folder_path.txt', 'w') as f:
        f.write(figure_folder_path)
    with open('model_folder_path.txt', 'w') as f:
        f.write(model_folder_path)
        
    # Print the paths for confirmation
    print(f"Figure folder path written to: figure_folder_path.txt")
    print(f"Model folder path written to: model_folder_path.txt")

