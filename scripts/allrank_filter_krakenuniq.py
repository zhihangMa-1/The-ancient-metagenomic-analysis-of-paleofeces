#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description: Filters KrakenUniq reports based on unique k-mer counts and taxonomic read depth.
Usage: python allrank_filter_krakenuniq.py <krakenuniq_report> <min_kmers> <min_tax_reads>
"""

import sys
import pandas as pd

def filter_krakenuniq():
    # 1. Validate Input Arguments
    if len(sys.argv) != 4:
        print("Usage: python allrank_filter_krakenuniq.py <krakenuniq_report> <min_kmers> <min_tax_reads>")
        sys.exit(1)

    krakenuniq_report = sys.argv[1]
    min_kmers = int(sys.argv[2])
    min_tax_reads = int(sys.argv[3])

    # 2. Load KrakenUniq Report
    # KrakenUniq outputs use tab-separated values. Comments starting with '#' are skipped.
    try:
        df = pd.read_csv(krakenuniq_report, sep='\t', comment='#')
    except Exception as e:
        print(f"Error: Failed to read the report file. {e}")
        sys.exit(1)

    # Clean column names (strip leading/trailing whitespace)
    df.columns = [col.strip() for col in df.columns]

    # Convert numeric columns to ensure proper filtering
    numeric_columns = ['%', 'reads', 'taxReads', 'kmers', 'dup', 'cov', 'taxID']
    for col in numeric_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    print(f"[*] Initial entries in report: {len(df)}")

    # 3. Apply Stringent Filters
    # Filter by Breadth of Coverage (Unique k-mers)
    df = df[df["kmers"] >= min_kmers]
    print(f"[*] Entries after k-mer filter (>= {min_kmers}): {len(df)}")

    # Filter by Depth of Coverage (Taxonomic Reads)
    df = df[df["taxReads"] >= min_tax_reads]
    print(f"[*] Entries after taxonomic read filter (>= {min_tax_reads}): {len(df)}")

    # 4. Split and Save by Taxonomic Rank
    # Standard taxonomic hierarchy used for downstream analysis
    ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

    for rank in ranks:
        rank_df = df[df["rank"] == rank]
        
        if not rank_df.empty:
            # Sort by relative abundance (%) in descending order
            rank_df = rank_df.sort_values(by="%", ascending=False)
            
            # Generate output filename (e.g., sample.report.species.filtered)
            output_path = f"{krakenuniq_report}.{rank}.filtered"
            rank_df.to_csv(output_path, sep="\t", index=False)
            print(f"[+] Saved filtered {rank} results to: {output_path}")

if __name__ == "__main__":
    filter_krakenuniq()
