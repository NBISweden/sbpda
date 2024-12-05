#!/usr/bin/env python
import argparse
import pandas as pd
import os

def validate_samplesheet(samplesheet_path, output_path):
    # Load samplesheet
    samplesheet = pd.read_csv(samplesheet_path)

    # Perform validation checks (e.g., file presence)
    samplesheet['adata_path'] = samplesheet['adata_path'].apply(lambda x: os.path.exists(x))

    # Save validated samplesheet
    samplesheet.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True, help="Path to samplesheet.csv")
    parser.add_argument("--output", default = 'validated_samplesheet.csv', required=False, help="Output path for validated samplesheet")
    args = parser.parse_args()

    validate_samplesheet(args.samplesheet, args.output)
