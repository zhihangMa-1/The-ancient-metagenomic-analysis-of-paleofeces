#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description: Enriches KrakenUniq filtered results with full taxonomic lineage (Superkingdom to Species).
Usage: python get_lineages_all.py <kraken_filtered_file> <output_file> <lineage_info_tsv>
"""

import pandas as pd
import sys
import os

def load_lineage_database(lineage_file):
    """Load the taxonomy database and create a lookup dictionary."""
    if not os.path.exists(lineage_file):
        print(f"Error: Taxonomy database not found at {lineage_file}")
        sys.exit(1)
        
    print(f"[*] Loading taxonomy database: {lineage_file}...")
    lineage_df = pd.read_csv(lineage_file, sep="\t")
    lineage_df.columns = lineage_df.columns.str.strip()
    
    # Using 'index' orientation to create a fast lookup dictionary: {taxID: {rank, name, parentTaxID}}
    return lineage_df.set_index('taxID').to_dict(orient='index')

def trace_full_lineage(tax_id, lineage_dict):
    """Recursively trace parent TaxIDs to construct the full lineage path."""
    path = {
        "superkingdom": None, "phylum": None, "class": None, 
        "order": None, "family": None, "genus": None, "species": None
    }
    
    current_id = tax_id
    # TaxID 1 is the root of the NCBI taxonomy tree
    while current_id in lineage_dict and current_id != 1:
        entry = lineage_dict[current_id]
        rank = entry["rank"]
        name = entry["name"]
        
        if rank in path:
            path[rank] = name
            
        current_id = entry["parentTaxID"]
        
    return path

def main():
    if len(sys.argv) < 3:
        print("Usage: python get_lineages_all.py <kraken_filtered_file> <output_file> [lineage_info_tsv]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Priority: 1. Argument from command line, 2. Default relative path, 3. Error
    if len(sys.argv) == 4:
        lineage_db_path = sys.argv[3]
    else:
        # Suggesting a standardized location for the database
        lineage_db_path = "./references/taxonomy/lineage_info.tsv"

    # 1. Load data
    print(f"[*] Reading filtered KrakenUniq results: {input_file}")
    kraken_df = pd.read_csv(input_file, sep="\t")
    
    if 'taxID' not in kraken_df.columns:
        print("Error: 'taxID' column not found in input file.")
        sys.exit(1)

    # 2. Load Lineage DB
    lineage_dict = load_lineage_database(lineage_db_path)

    # 3. Process each TaxID
    print("[*] Tracing lineages for all identified taxa...")
    enriched_results = []
    for _, row in kraken_df.iterrows():
        tax_id = row['taxID']
        lineage_path = trace_full_lineage(tax_id, lineage_dict)
        # Merge the original row data with the new lineage columns
        enriched_row = {**row, **lineage_path}
        enriched_results.append(enriched_row)

    # 4. Save results
    final_df = pd.DataFrame(enriched_results)
    final_df.to_csv(output_file, sep="\t", index=False)
    print(f"[+] Enrichment complete. Results saved to: {output_file}")

if __name__ == "__main__":
    main()
