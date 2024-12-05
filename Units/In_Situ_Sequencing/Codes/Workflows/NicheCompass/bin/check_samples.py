#!/usr/bin/env python
import argparse
import pandas as pd
import os
import sys

def validate_samplesheet(samplesheet_path, data_path, output_path):
    # Load samplesheet
    samplesheet = pd.read_csv(samplesheet_path)

    # Get the absolute path of the data directory
    data_dir = os.path.abspath(data_path)

    # Perform validation checks (e.g., file presence)
    samplesheet['adata_path'] = samplesheet['adata_path'].apply(
        lambda x: os.path.join(data_dir, x) if os.path.exists(os.path.join(data_dir, x)) else False
    )

    # Check if any files exist
    if not any(samplesheet['adata_path']):
        print("Error: None of the specified files were found.")
        sys.exit(1)

    # Remove rows where files don't exist
    valid_entries = samplesheet['adata_path'].astype(bool).sum()
    samplesheet = samplesheet[samplesheet['adata_path'].astype(bool)]

    # Save validated samplesheet
    samplesheet.to_csv(output_path, index=False)
    print(f"Validated samplesheet saved to {output_path}")
    print(f"Found {valid_entries} valid entries out of {len(samplesheet)} total.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True, help="Path to samplesheet.csv")
    parser.add_argument("--data_path", required=True, help="Path to the directory containing the data files")
    parser.add_argument("--output", default='validated_samplesheet.csv', required=False, help="Output path for validated samplesheet")
    args = parser.parse_args()

    validate_samplesheet(args.samplesheet, args.data_path, args.output)