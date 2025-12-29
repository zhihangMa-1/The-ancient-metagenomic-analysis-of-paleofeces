#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description: Consolidates individual sample taxonomic reports into a single abundance matrix.
Usage: python generate_abundance_matrix.py <input_dir> <sample_list_file> <output_csv>
"""

import pandas as pd
import glob
import os
import sys

def main():
    # 1. Parameter Validation
    if len(sys.argv) != 4:
        print("Usage: python generate_abundance_matrix.py <input_dir> <sample_list_file> <output_csv>")
        sys.exit(1)

    input_dir = sys.argv[1]
    sample_order_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Expected filename suffix from the lineage enrichment step
    FILE_SUFFIX = 'krakenuniq.output.species.filtered.with_lineage.tsv'

    # 2. Load Targeted Sample Order
    try:
        with open(sample_order_file, 'r') as f:
            sample_order = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error: Failed to read sample list '{sample_order_file}': {e}")
        sys.exit(1)

    # 3. Locate Result Files
    # Search recursively for enriched taxonomic files
    search_pattern = os.path.join(input_dir, '*/' + FILE_SUFFIX)
    tsv_files = glob.glob(search_pattern)

    if not tsv_files:
        print(f"Error: No files matching '{FILE_SUFFIX}' found in '{input_dir}'.")
        sys.exit(1)

    print(f"[*] Found {len(tsv_files)} potential result files. Processing...")

    # 4. Aggregate Data
    combined_list = []
    for file in tsv_files:
        # The sample name is derived from the parent directory name
        sample_name = os.path.basename(os.path.dirname(file))
        
        # Only include samples specified in the order file
        if sample_name not in sample_order:
            continue
            
        try:
            df = pd.read_csv(file, sep='\t')
            df['Sample'] = sample_name
            combined_list.append(df)
        except Exception as e:
            print(f"Warning: Skipping file '{file}' due to error: {e}")
            continue

    if not combined_list:
        print("Error: No data matched the provided sample list.")
        sys.exit(1)

    combined_df = pd.concat(combined_list, ignore_index=True)

    # 5. Define Metadata and Reshape Matrix
    # Columns to be used as taxonomic identifiers
    taxa_metadata = [
        'taxID', 'superkingdom', 'phylum', 'class', 
        'order', 'family', 'genus', 'species'
    ]
    
    # Verify presence of required columns
    required_cols = taxa_metadata + ['Sample', 'reads']
    missing = [c for c in required_cols if c not in combined_df.columns]
    if missing:
        print(f"Error: Missing required columns in input files: {missing}")
        sys.exit(1)

    # Pivot the data: Rows = Taxa, Columns = Samples, Values = Reads
    abundance_matrix = combined_df.pivot_table(
        index=taxa_metadata,
        columns='Sample',
        values='reads',
        fill_value=0
    )

    # 6. Reorder and Export
    # Match the column order to the provided sample list
    available_samples = [s for s in sample_order if s in abundance_matrix.columns]
    missing_samples = [s for s in sample_order if s not in abundance_matrix.columns]
    
    if missing_samples:
        print(f"[*] Note: The following samples were listed but not found in data: {missing_samples}")

    abundance_matrix = abundance_matrix[available_samples]

    # Save to CSV (compatible with R/Phyloseq)
    abundance_matrix.to_csv(output_file)
    print(f"[+] Success! Abundance matrix saved to '{output_file}'")
    print(f"[*] Final matrix dimensions: {abundance_matrix.shape}")

if __name__ == "__main__":
    main()
