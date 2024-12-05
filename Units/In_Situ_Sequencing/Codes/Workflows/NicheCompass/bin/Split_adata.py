#!/usr/bin/env python
import argparse
import os
import scanpy as sc

def split_adata_by_column(input_file, column_name, output_dir):
    """
    Splits a merged AnnData file into individual files based on a specified column.

    Parameters:
    input_file (str): Path to the merged AnnData file.
    column_name (str): Column in adata.obs to split the data by.
    output_dir (str): Directory to save the individual files.

    Returns:
    None
    """
    # Load the merged AnnData file
    print(f"Loading AnnData file: {input_file}")
    combined_adata = sc.read_h5ad(input_file)
    
    if column_name not in combined_adata.obs.columns:
        print(f"Error: Column '{column_name}' not found in adata.obs.")
        return

    # Get unique values in the specified column
    unique_values = combined_adata.obs[column_name].unique()
    print(f"Unique values in column '{column_name}': {unique_values}")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop through each unique value and save a subset of the data
    for value in unique_values:
        print(f"Processing subset: {value}")
        
        # Subset the data for the current value
        adata_individual = combined_adata[combined_adata.obs[column_name] == value].copy()
        
        # Define the output file name
        output_file = os.path.join(output_dir, f"adata_{column_name}_{value}.h5ad")
        
        # Save the individual AnnData object
        adata_individual.write(output_file)
        print(f"Saved: {output_file}")

    print("Splitting complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a merged AnnData file into individual files based on a specified column.")
    parser.add_argument('--input_file', required=True, help="Path to the merged AnnData file.")
    parser.add_argument('--column_name', required=True, help="Column in adata.obs to split the data by.")
    parser.add_argument('--output_dir', required=True, help="Directory to save the individual files.")
    
    args = parser.parse_args()

    split_adata_by_column(args.input_file, args.column_name, args.output_dir)
