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

def prepare_folders(base_path, data_path, tempDir, species, samplesheet):
    if not samplesheet or not os.path.exists(samplesheet):
        print("Error: samplesheet is missing or does not exist. Please provide a valid CSV file with path information.")
        sys.exit(1)
    # Check and create base, data, and temporary directories if they do not exist
    for path, label in zip([base_path, data_path, tempDir], ['Base directory', 'Data directory', 'Temporary directory']):
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"{label} created at: {path}")
        else:
            print(f"{label} exists at: {path}")

    
    # Create folders for gene programs and gene annotations
    ga_data_folder_path = os.path.join(data_path, 'gene_annotations')
    gp_data_folder_path = os.path.join(data_path, 'gene_programs')
    artifacts_folder_path = os.path.join(tempDir, 'artifacts', 'sample_integration')

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
    model_folder_path = os.path.join(artifacts_folder_path, 'model')
    figure_folder_path = os.path.join(artifacts_folder_path, 'figures')
    os.makedirs(model_folder_path, exist_ok=True)
    os.makedirs(figure_folder_path, exist_ok=True)

    # Print paths for confirmation
    print(f"Gene annotation folder: {ga_data_folder_path}")
    print(f"Gene programs folder: {gp_data_folder_path}")
    print(f"Artifacts folder: {artifacts_folder_path}")
    print(f"Model folder: {model_folder_path}")
    print(f"Figures folder: {figure_folder_path}")

    return('figure_folder_path' : figure_folder_path, 
           'model_folder_path'  : model_folder_path,)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare necessary folders for NicheCompass analysis.")
    parser.add_argument('--base_path', required=True, help="Base directory for the project")
    parser.add_argument('--data_path', required=True, help="Directory for storing gene annotations and programs data")
    parser.add_argument('--tempDir', required=True, help="Temporary directory for storing intermediate files")
    parser.add_argument('--species', default="human", help="Species for NicheCompass analysis (default: human)")
    parser.add_argument('--samplesheet', help="Sample sheet file")
    args = parser.parse_args()

    # Run folder preparation and retrieve paths
    figure_folder_path, model_folder_path = prepare_folders(
        base_path=args.base_path,
        data_path=args.data_path,
        tempDir=args.tempDir,
        species=args.species,
        samplesheet=args.samplesheet
    )

    # Print the paths for confirmation
    print(f"Returned Figure Folder Path: {figure_folder_path}")
    print(f"Returned Model Folder Path: {model_folder_path}")