#!/usr/bin/env python
import argparse
import os
import random
import warnings
from datetime import datetime
import time

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

# gp extraction
lrt_interactions_version = "v2"
keep_target_genes_ratio_val=1.0
max_n_target_genes_per_gp_val=250
min_genes_per_gp_val=2
min_source_genes_per_gp_val=1
min_target_genes_per_gp_val=1
max_genes_per_gp_val=None
max_source_genes_per_gp_val=None
max_target_genes_per_gp_val=None

# Architecture
cat_covariates_embeds_injection = ["gene_expr_decoder"]
cat_covariates_no_edges = [True]

# Analysis 
differential_gp_test_results_key = 'nichecompass_differential_gp_test_results'


def create_checkpoint(file_path):
    """Create a checkpoint file to indicate the process has completed."""
    with open(file_path, 'w') as f:
        f.write("Completed")

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

    return figure_folder_path, model_folder_path, omnipath_lr_network_file_path, \
        nichenet_lr_network_file_path, nichenet_ligand_target_matrix_file_path, \
        gene_orthologs_mapping_file_path, mebocost_enzyme_sensor_interactions_folder_path

def load_anndata(samplesheet_path, data_path, spatial_key, n_neighbors, adj_key, model_folder_path):
    valid_files_cntr = 0
    missing_files_cntr = 0
    missing_files = []
    valid_files = []
    adata_batch_list = []
    
    samples = pd.read_csv(samplesheet_path)
    for fl in range(samples.shape[0]):
        file_path = samples['adata_path'][fl]
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
        # Check if 'spatial' exists in obsm
        if 'spatial' not in adata_batch.obsm:
            print("'spatial' key missing in adata.obsm. Attempting to create it...")

            # Check for the necessary columns in adata.obs
            coord_keys = ('x_centroid', 'y_centroid', 'spatial')
            if coord_keys[0] in adata_batch.obs and coord_keys[1] in adata_batch.obs:
                adata_batch.obsm[coord_keys[2]] = pd.concat(
                    [adata_batch.obs[coord_keys[0]], adata_batch.obs[coord_keys[1]]], axis=1
                ).to_numpy()
                print(f"'spatial' key successfully created in adata.obsm.")
            else:
                raise ValueError(
                    f"Cannot create 'spatial' key: Missing '{coord_keys[0]}' or '{coord_keys[1]}' in adata.obs."
                )

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
    return adata, adata_batch_list

def combine_spatial_neighborhood_graphs(adata, adata_batch_list):
    """
    Combine spatial neighborhood graphs as disconnected components.
    
    Parameters:
    adata (AnnData): The combined AnnData object.
    adata_batch_list (list): List of AnnData objects.
    adj_key has default value set in the global variable.

    Returns:
    AnnData: The input AnnData object with updated obsp[adj_key].
    """
    batch_connectivities = []
    len_before_batch = 0
    
    for i in range(len(adata_batch_list)):
        if i == 0:  # first batch
            after_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[0].shape[0],
                (adata.shape[0] - adata_batch_list[0].shape[0])))
            batch_connectivities.append(sp.hstack(
                (adata_batch_list[0].obsp[adj_key],
                after_batch_connectivities_extension)))
        elif i == (len(adata_batch_list) - 1):  # last batch
            before_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0],
                (adata.shape[0] - adata_batch_list[i].shape[0])))
            batch_connectivities.append(sp.hstack(
                (before_batch_connectivities_extension,
                adata_batch_list[i].obsp[adj_key])))
        else:  # middle batches
            before_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0], len_before_batch))
            after_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0],
                (adata.shape[0] - adata_batch_list[i].shape[0] - len_before_batch)))
            batch_connectivities.append(sp.hstack(
                (before_batch_connectivities_extension,
                adata_batch_list[i].obsp[adj_key],
                after_batch_connectivities_extension)))
        
        len_before_batch += adata_batch_list[i].shape[0]
    
    adata.obsp[adj_key] = sp.vstack(batch_connectivities)
    sc.write(os.path.join(model_folder_path, "adata_merged_spatial_neighbour_graph.h5ad"), adata)
    return adata


def process_gene_programs(
    species,
    figure_folder_path,
    omnipath_lr_network_file_path,
    nichenet_lr_network_file_path,
    nichenet_ligand_target_matrix_file_path,
    gene_orthologs_mapping_file_path,
    mebocost_enzyme_sensor_interactions_folder_path,
    adata,
    gp_targets_mask_key,
    gp_targets_categories_mask_key,
    gp_sources_mask_key,
    gp_sources_categories_mask_key,
    gp_names_key,
):
    """Processes and combines gene programs from OmniPath, NicheNet, and MEBOCOST."""
    # Retrieve OmniPath GPs
    print("Downloading OmniPath gene programs....")
    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species=species,
        load_from_disk=False,
        save_to_disk=True,
        lr_network_file_path=omnipath_lr_network_file_path,
        gene_orthologs_mapping_file_path="",
        plot_gp_gene_count_distributions=False,
        gp_gene_count_distributions_save_path=f"{figure_folder_path}/omnipath_gp_gene_count_distributions.svg",
    )

    # Retrieve NicheNet GPs
    print("Downloading NicheNet gene programs....")
    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species=species,
        version=lrt_interactions_version,
        keep_target_genes_ratio=keep_target_genes_ratio_val,
        max_n_target_genes_per_gp=max_n_target_genes_per_gp_val,
        load_from_disk=False,
        save_to_disk=True,
        lr_network_file_path=nichenet_lr_network_file_path,
        ligand_target_matrix_file_path=nichenet_ligand_target_matrix_file_path,
        gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
        plot_gp_gene_count_distributions=False,
    )

    # Retrieve MEBOCOST GPs
    print("Downloading Mebocost gene programs....")
    mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
        dir_path=mebocost_enzyme_sensor_interactions_folder_path,
        species=species,
        plot_gp_gene_count_distributions=False,
    )

    # Filter and combine GPs
    print("Filtering and combining GPs...")
    gp_dicts = [omnipath_gp_dict, nichenet_gp_dict, mebocost_gp_dict]
    combined_gp_dict = filter_and_combine_gp_dict_gps_v2(gp_dicts, verbose=False)

    # Add the GP dictionary as binary masks to the adata
    print("Adding GP dictionary to AnnData object...")
    add_gps_from_gp_dict_to_adata(
        gp_dict=combined_gp_dict,
        adata=adata,
        gp_targets_mask_key=gp_targets_mask_key,
        gp_targets_categories_mask_key=gp_targets_categories_mask_key,
        gp_sources_mask_key=gp_sources_mask_key,
        gp_sources_categories_mask_key=gp_sources_categories_mask_key,
        gp_names_key=gp_names_key,
        min_genes_per_gp=min_genes_per_gp_val,
        min_source_genes_per_gp=min_source_genes_per_gp_val,
        min_target_genes_per_gp=min_target_genes_per_gp_val,
        max_genes_per_gp=max_genes_per_gp_val,
        max_source_genes_per_gp=max_source_genes_per_gp_val,
        max_target_genes_per_gp=max_target_genes_per_gp_val,
    )
    sc.write(os.path.join(model_folder_path, "adata_merged_spatial_neighbour_graph_gps.h5ad"), adata)
    return adata

def initialize_model(adata, use_cuda_if_available, adj_key, 
                     cat_covariates_embeds_injection, cat_covariates_keys, 
                     cat_covariates_no_edges, cat_covariates_embeds_nums, 
                     gp_names_key, active_gp_names_key, gp_targets_mask_key, 
                     gp_targets_categories_mask_key, gp_sources_mask_key, 
                     gp_sources_categories_mask_key, latent_key, 
                     conv_layer_encoder, active_gp_thresh_ratio):
    """
    Initializes the NicheCompass model.

    Parameters:
    adata (AnnData): The AnnData object.
    use_cuda_if_available (bool): Whether to use CUDA if available.
    Other parameters: Model initialization configurations.

    Returns:
    model: The initialized NicheCompass model.
    """
    # Convert to list if passed as a single string
    if isinstance(cat_covariates_keys, str):
        cat_covariates_keys = [cat_covariates_keys]
    if isinstance(cat_covariates_embeds_nums, str):
        cat_covariates_embeds_nums = [int(cat_covariates_embeds_nums)]

    print('Initializing the model....')
    adata.X = sp.csr_matrix(adata.X)

    # Device setup
    device = torch.device("cuda" if use_cuda_if_available else "cpu")
    print(device)

    # Initialize the model
    model = NicheCompass(
        adata=adata,
        counts_key=None,
        adj_key=adj_key,
        cat_covariates_embeds_injection=cat_covariates_embeds_injection,
        cat_covariates_keys=cat_covariates_keys,
        cat_covariates_no_edges=cat_covariates_no_edges,
        cat_covariates_embeds_nums=cat_covariates_embeds_nums,
        gp_names_key=gp_names_key,
        active_gp_names_key=active_gp_names_key,
        gp_targets_mask_key=gp_targets_mask_key,
        gp_targets_categories_mask_key=gp_targets_categories_mask_key,
        gp_sources_mask_key=gp_sources_mask_key,
        gp_sources_categories_mask_key=gp_sources_categories_mask_key,
        latent_key=latent_key,
        conv_layer_encoder=conv_layer_encoder,
        active_gp_thresh_ratio=active_gp_thresh_ratio,
    )  

    return model

def train_model_and_identify_niches(model, latent_key, model_folder_path, n_epochs, 
                                    n_epochs_all_gps, lr, lambda_edge_recon, 
                                    lambda_gene_expr_recon, lambda_l1_masked, 
                                    edge_batch_size, n_sampled_neighbors, 
                                    use_cuda_if_available, latent_leiden_resolution, 
                                    latent_cluster_key):
    """
    Trains the NicheCompass model, computes latent clusters, and identifies niches.

    Parameters:
    model: The NicheCompass model to be trained.
    latent_key (str): Key for latent representation in AnnData.
    model_folder_path (str): Path to save the trained model.
    latent_leiden_resolution (float): Resolution parameter for Leiden clustering.
    latent_cluster_key (str): Key to store the clustering results in adata.obs.
    Other parameters: Training configurations.

    Returns:
    tuple: Updated model and a dictionary mapping cluster labels to colors.
    """
    print('Training the model....')

    # Train the model
    model.train(
        n_epochs=n_epochs,
        n_epochs_all_gps=n_epochs_all_gps,
        lr=lr,
        lambda_edge_recon=lambda_edge_recon,
        lambda_gene_expr_recon=lambda_gene_expr_recon,
        lambda_l1_masked=lambda_l1_masked,
        edge_batch_size=edge_batch_size,
        n_sampled_neighbors=n_sampled_neighbors,
        use_cuda_if_available=use_cuda_if_available,
        verbose=False,
    )
    print("Model training completed.")

    # Compute latent neighbor graph
    print("Computing latent neighbor graph...")
    sc.pp.neighbors(model.adata, use_rep=latent_key, key_added=latent_key)

    # Compute UMAP embedding
    print("Computing UMAP embedding...")
    sc.tl.umap(model.adata, neighbors_key=latent_key)

    # Perform Leiden clustering
    print(f"Performing Leiden clustering with resolution={latent_leiden_resolution}...")
    sc.tl.leiden(
        adata=model.adata,
        resolution=latent_leiden_resolution,
        key_added=latent_cluster_key,
        neighbors_key=latent_key
    )
    print(f"Leiden clustering completed. Results stored in '{latent_cluster_key}'.")

    # Create a color dictionary for the clusters
    latent_cluster_colors = create_new_color_dict(
        adata=model.adata,
        cat_key=latent_cluster_key
    )
    print(f"Color dictionary for latent clusters created.")

    # Save the trained model and updated AnnData
    model.save(dir_path=model_folder_path, overwrite=True, save_adata=True, adata_file_name="adata.h5ad")
    print(f"Model and AnnData saved in {model_folder_path}")

    return model, latent_cluster_colors


def plot_batches_latent_physical_space(model, sample_key, cat_covariates_keys, 
                                       batch_colors, figure_folder_path, 
                                       spot_size, save_fig=True):
    """
    Plots batches in latent and physical space.


    Parameters:
    model: The trained NicheCompass model.
    sample_key (str): Key for sample information in AnnData.
    cat_covariates_keys (list): List of categorical covariates keys.
    batch_colors (dict): Dictionary of batch colors.
    figure_folder_path (str): Path to save the figure.
    spot_size (int): Spot size for spatial plots.
    save_fig (bool): Whether to save the figure.

    Returns:
    None
    """
    samples = model.adata.obs[sample_key].unique().tolist()
    file_path = f"{figure_folder_path}/batches_latent_physical_space.png"
    fig = plt.figure(figsize=(12, 14))
    title = fig.suptitle(t="NicheCompass Batches in Latent and Physical Space",
                         y=0.96, x=0.55, fontsize=20)
    spec1 = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[3, 2])
    spec2 = gridspec.GridSpec(ncols=len(samples), nrows=2, height_ratios=[3, 2])
    axs = []
    axs.append(fig.add_subplot(spec1[0]))
    sc.pl.umap(adata=model.adata,
               color=[cat_covariates_keys[0]],
               palette=batch_colors,
               title="Batches in Latent Space",
               ax=axs[0],
               show=False)
    for idx, sample in enumerate(samples):
        axs.append(fig.add_subplot(spec2[len(samples) + idx]))
        sc.pl.spatial(adata=model.adata[model.adata.obs[sample_key] == sample],
                      color=[cat_covariates_keys[0]],
                      palette=batch_colors,
                      spot_size=spot_size,
                      title=f"Batches in Physical Space \n(Sample: {sample})",
                      ax=axs[idx+1],
                      show=False)

    # Create shared legend
    handles, labels = axs[0].get_legend_handles_labels()
    lgd = fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.98, 0.5))
    axs[0].get_legend().remove()

    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    if save_fig:
        fig.savefig(file_path, bbox_extra_artists=(lgd, title), bbox_inches="tight")

def plot_niches_latent_physical_space(model, sample_key, latent_cluster_key, 
                                      latent_cluster_colors, figure_folder_path, 
                                      spot_size, latent_leiden_resolution, save_fig=True):
    
    """
    Plots niches in latent and physical space.

    Parameters:
    model: The trained NicheCompass model.
    sample_key (str): Key for sample information in AnnData.
    latent_cluster_key (str): Key for latent cluster information in AnnData.
    latent_cluster_colors (dict): Dictionary of latent cluster colors.
    figure_folder_path (str): Path to save the figure.
    spot_size (int): Spot size for spatial plots.
    latent_leiden_resolution (float): Leiden resolution for latent clusters.
    save_fig (bool): Whether to save the figure.

    Returns:
    None
    """ 
    samples = model.adata.obs[sample_key].unique().tolist()
    file_path = f"{figure_folder_path}/res_{latent_leiden_resolution}_niches_latent_physical_space.png"
    fig = plt.figure(figsize=(12, 14))
    title = fig.suptitle(t="NicheCompass Niches in Latent and Physical Space",
                         y=0.96, x=0.55, fontsize=20)
    spec1 = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[3, 2])
    spec2 = gridspec.GridSpec(ncols=len(samples), nrows=2, height_ratios=[3, 2])
    axs = []
    axs.append(fig.add_subplot(spec1[0]))
    sc.pl.umap(adata=model.adata,
               color=[latent_cluster_key],
               palette=latent_cluster_colors,
               title="Niches in Latent Space",
               ax=axs[0],
               show=False)
    for idx, sample in enumerate(samples):
        axs.append(fig.add_subplot(spec2[len(samples) + idx]))
        sc.pl.spatial(adata=model.adata[model.adata.obs[sample_key] == sample],
                      color=[latent_cluster_key],
                      palette=latent_cluster_colors,
                      spot_size=spot_size,
                      title=f"(Sample: {sample})",
                      ax=axs[idx+1],
                      show=False)

    # Create shared legend
    handles, labels = axs[0].get_legend_handles_labels()
    lgd = fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.98, 0.5))
    axs[0].get_legend().remove()

    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    if save_fig:
        fig.savefig(file_path, bbox_extra_artists=(lgd, title), bbox_inches="tight")

def plot_niche_composition(model, latent_cluster_key, cell_type_key, 
                           figure_folder_path, latent_leiden_resolution, save_fig=True):
    """
    Plots the cell type composition of niches.
    
    Parameters:
    model: The trained NicheCompass model.
    latent_cluster_key (str): Key for latent cluster information in AnnData.
    cell_type_key (str): Key for cell type information in AnnData.
    figure_folder_path (str): Path to save the figure.
    latent_leiden_resolution (float): Leiden resolution for latent clusters.
    save_fig (bool): Whether to save the figure.
    
    Returns:
    None
    """

    file_path = f"{figure_folder_path}/res_{latent_leiden_resolution}_niche_composition_cell_types.png"
    df_counts = (model.adata.obs.groupby([latent_cluster_key, cell_type_key])
                 .size().unstack())
    df_counts.plot(kind="bar", stacked=True, figsize=(10, 10))
    legend = plt.legend(bbox_to_anchor=(1, 1), loc="upper left", prop={'size': 10})
    legend.set_title("Cell Type Annotations", prop={'size': 10})
    plt.title("Cell Type Composition of Niches")
    plt.xlabel("Niche")
    plt.ylabel("Cell Counts")
    if save_fig:
        plt.savefig(file_path, bbox_extra_artists=(legend,), bbox_inches="tight")


def add_groupid_for_tissuumaps(adata, group_by_key="replicate", groupid_column="groupid", enable_tissuumaps=False):
    """
    Adds a 'groupid' column to AnnData for TissUUmaps visualization.
    
    Parameters:
    adata (AnnData): The AnnData object to modify.
    group_by_key (str): The column in adata.obs to group by (default: "replicate").
    groupid_column (str): The name of the new column to add (default: "groupid").
    enable_tissuumaps (bool): Whether to enable this visualization preparation step.
    
    Returns:
    adata (AnnData): Modified AnnData object with the new 'groupid' column.
    """
    if enable_tissuumaps:
        print(f"Adding '{groupid_column}' column for TissUUmaps visualization based on '{group_by_key}'...")

        # Get unique values and sort them alphabetically
        unique_samples = sorted(adata.obs[group_by_key].unique())
        print(f"Unique values in '{group_by_key}': {unique_samples}")

        # Create mapping from unique values to integers
        sample_mapping = {sample: i for i, sample in enumerate(unique_samples)}
        print(f"Mapping: {sample_mapping}")

        # Add the groupid column
        adata.obs[groupid_column] = adata.obs[group_by_key].map(sample_mapping)

        print(f"'{groupid_column}' column added to adata.obs.")
    else:
        print("TissUUmaps visualization preparation is disabled. Skipping groupid creation.")
    
    return adata


# Separate the data preparation and analysis steps:
    # [x] samplesheet
    # [x] folder preparation
    # [x] loading the anndata objects
    # [x] Combine spatial neighborhood graphs as disconnected components
    # [x] Extract gene programs from gene-gene interaction networks
    # [x] Filter and combine gene programs
        # **** metabolite_enzyme_sensor_gps are delivered with NicheCompass package which means that these data which are tsv files for both human and mouse should be copied to new local data directories in each run. For this we can 
    # [x] Train the NicheCompass model
    # [x] Analyze the trained model
    # [x] Visualize the results
    # [x] Save the results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare necessary folders for NicheCompass analysis.")
    parser.add_argument('--base_path', required=True, help="Base directory for the project")
    parser.add_argument('--data_path', required=True, help="Directory for storing gene annotations and programs data")
    parser.add_argument('--tempDir', required=True, help="Temporary directory for storing intermediate files")
    parser.add_argument('--species', default="human", help="Species for NicheCompass analysis (default: human)")
    parser.add_argument('--samplesheet', help="Sample sheet file")

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
    parser.add_argument('--sample_key', help="Sample key")
    parser.add_argument('--spot_size', type=int, help="Spot size")

    # TissUUmaps
    parser.add_argument('--enable_tissuumaps', action='store_true',help="Enable preparation of data for TissUUmaps visualization by adding a 'groupid' column.")

    # Additional parameters
    parser.add_argument('--spatial_key', help="Spatial key for analysis")
    parser.add_argument('--n_neighbors', type=int, help="Number of neighbors")
    parser.add_argument('--mebocost_dir', default="../assets/gene_programs/metabolite_enzyme_sensor_gps/", help="Directory containing MEBOCOST gene program files")
    
    args = parser.parse_args()


    start_time = time.time()
    
    # Convert to list if passed as a single string
    if isinstance(args.cat_covariates_keys, str):
        args.cat_covariates_keys = [args.cat_covariates_keys]
        
    # Run folder preparation and retrieve paths
    figure_folder_path, model_folder_path, omnipath_lr_network_file_path, \
        nichenet_lr_network_file_path, nichenet_ligand_target_matrix_file_path, \
        gene_orthologs_mapping_file_path, mebocost_enzyme_sensor_interactions_folder_path = prepare_folders(
        base_path=args.base_path,
        data_path=args.data_path,
        tempDir=args.tempDir,
        species=args.species
    )

    samplesheet_path = load_samplesheet(args.samplesheet)
    print(f'Loding anndata objects from {samplesheet_path}...')
    adata, adata_batch_list = load_anndata(
        samplesheet_path=args.samplesheet,
        data_path=args.data_path,
        spatial_key=args.spatial_key,
        n_neighbors=args.n_neighbors,
        adj_key=adj_key,
        model_folder_path = model_folder_path)
    print(f"adata is saved in {model_folder_path}\n")

    print(f"Combined spatial neighborhood graphs as disconnected components...")
    adata = combine_spatial_neighborhood_graphs(adata, adata_batch_list)
    
    print(f"Processing gene programs...")
    adata = process_gene_programs(
        species=args.species,
        figure_folder_path=figure_folder_path,
        omnipath_lr_network_file_path=omnipath_lr_network_file_path,
        nichenet_lr_network_file_path=nichenet_lr_network_file_path,
        nichenet_ligand_target_matrix_file_path=nichenet_ligand_target_matrix_file_path,
        gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
        mebocost_enzyme_sensor_interactions_folder_path=args.mebocost_dir,
        adata=adata,
        gp_targets_mask_key=gp_targets_mask_key,
        gp_targets_categories_mask_key=gp_targets_categories_mask_key,
        gp_sources_mask_key=gp_sources_mask_key,
        gp_sources_categories_mask_key=gp_sources_categories_mask_key,
        gp_names_key=gp_names_key,
    )


    # Initialize the model
    model = initialize_model(
        adata=adata,
        use_cuda_if_available=args.use_cuda_if_available,
        adj_key=adj_key,
        cat_covariates_embeds_injection=cat_covariates_embeds_injection,
        cat_covariates_keys=args.cat_covariates_keys,
        cat_covariates_no_edges=cat_covariates_no_edges,
        cat_covariates_embeds_nums=args.cat_covariates_embeds_nums,
        gp_names_key=gp_names_key,
        active_gp_names_key=active_gp_names_key,
        gp_targets_mask_key=gp_targets_mask_key,
        gp_targets_categories_mask_key=gp_targets_categories_mask_key,
        gp_sources_mask_key=gp_sources_mask_key,
        gp_sources_categories_mask_key=gp_sources_categories_mask_key,
        latent_key=latent_key,
        conv_layer_encoder=args.conv_layer_encoder,
        active_gp_thresh_ratio=args.active_gp_thresh_ratio,
    )

    
    # Train the model
    latent_cluster_key = f"latent_leiden_{str(args.latent_leiden_resolution)}"
    model, latent_cluster_colors = train_model_and_identify_niches(
        model=model,
        latent_key=latent_key,
        model_folder_path=model_folder_path,
        n_epochs=args.n_epochs,
        n_epochs_all_gps=args.n_epochs_all_gps,
        lr=args.lr,
        lambda_edge_recon=args.lambda_edge_recon,
        lambda_gene_expr_recon=args.lambda_gene_expr_recon,
        lambda_l1_masked=args.lambda_l1_masked,
        edge_batch_size=args.edge_batch_size,
        n_sampled_neighbors=args.n_sampled_neighbors,
        use_cuda_if_available=args.use_cuda_if_available,
        latent_leiden_resolution=args.latent_leiden_resolution,
        latent_cluster_key=latent_cluster_key
    )
    
    # Plot batches in latent and physical space
    batch_colors = create_new_color_dict(adata=model.adata, cat_key=args.cat_covariates_keys[0])
    plot_batches_latent_physical_space(
        model=model,
        sample_key=args.sample_key,
        cat_covariates_keys=args.cat_covariates_keys,
        batch_colors=batch_colors,
        figure_folder_path=figure_folder_path,
        spot_size=args.spot_size
    )

    # Plot niches in latent and physical space
    plot_niches_latent_physical_space(
        model=model,
        sample_key=args.sample_key,
        latent_cluster_key=latent_cluster_key,
        latent_cluster_colors=latent_cluster_colors,
        figure_folder_path=figure_folder_path,
        spot_size=args.spot_size,
        latent_leiden_resolution=args.latent_leiden_resolution
    )

    # Plot niche composition
    plot_niche_composition(
        model=model,
        latent_cluster_key=latent_cluster_key,
        cell_type_key=args.cell_type_key,
        figure_folder_path=figure_folder_path,
        latent_leiden_resolution=args.latent_leiden_resolution
    )

    # TissUUmaps
    if args.enable_tissuemaps:
        model.adata = add_groupid_for_tissuumaps(
        adata=model.adata,
        group_by_key="replicate",       # Column to group by
        groupid_column="groupid",       # New column name
        enable_tissuumaps=True  # User-specified flag
)

    
    # Save the trained model
    model.save(dir_path=model_folder_path, overwrite=True, save_adata=True, adata_file_name="adata.h5ad")
    print(f"Model saved in {model_folder_path}")

    
    # Finishing the analysis
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Analysis completed in {int(hours):02}:{int(minutes):02}:{int(seconds):02} (hh:mm:ss)")
