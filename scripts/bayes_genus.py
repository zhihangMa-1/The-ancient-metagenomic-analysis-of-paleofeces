#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description: Calculates Bayesian posterior probabilities for Genus-level assignments 
             based on BLAST bitscores, alignment length, and identity.
Usage: python bayes_genus_refinement.py <blast.tsv> <contig_lengths.tsv> <nodes.dmp> <names.dmp> <posterior_out.tsv> <abundance_out.tsv>
"""

import sys
from collections import defaultdict
import os

def main():
    if len(sys.argv) != 7:
        print("Usage: python bayes_genus_refinement.py <blast.tsv> <contig_lengths.tsv> <nodes.dmp> <names.dmp> <posterior_out.tsv> <abundance_out.tsv>")
        sys.exit(1)

    blast_file = sys.argv[1]
    length_file = sys.argv[2]
    nodes_file = sys.argv[3]
    names_file = sys.argv[4]
    posterior_file = sys.argv[5]
    abundance_file = sys.argv[6]

    # --- 1. Load Contig Lengths ---
    contig_len = {}
    with open(length_file) as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) >= 2:
                contig_len[cols[0]] = int(cols[1])

    # --- 2. Load NCBI Taxonomy Database ---
    parent = {}
    rank = {}
    print("[*] Parsing NCBI nodes.dmp...")
    with open(nodes_file) as f:
        for line in f:
            cols = line.strip().split("\t|\t")
            tid, parent_tid, r = cols[0].strip(), cols[1].strip(), cols[2].strip()
            parent[tid] = parent_tid
            rank[tid] = r

    taxid2name = {}
    print("[*] Parsing NCBI names.dmp...")
    with open(names_file) as f:
        for line in f:
            if "scientific name" not in line:
                continue
            cols = line.strip().split("\t|\t")
            tid = cols[0].strip()
            name = cols[1].strip()
            taxid2name[tid] = name

    def get_rank_tid(tid, target_rank):
        seen = set()
        while tid != "1" and tid in parent and tid not in seen:
            seen.add(tid)
            if rank.get(tid) == target_rank:
                return tid
            tid = parent[tid]
        return None

    # --- 3. Process BLAST Hits and Weight Scores ---
    # Score = Bitscore * (Alignment_Len / Contig_Len) * (Identity / 100)
    contig_genus_scores = defaultdict(lambda: defaultdict(float))
    contig_best_genus_tid = {}

    print(f"[*] Analyzing BLAST results: {blast_file}...")
    with open(blast_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 6: continue
            
            contig, tid = cols[0], cols[1]
            try:
                bitscore, align_len, pident = float(cols[2]), float(cols[3]), float(cols[4])
            except ValueError: continue

            if contig not in contig_len or contig_len[contig] == 0: continue
            
            # Weighted score calculation
            fraction = min(align_len / contig_len[contig], 1.0)
            score = bitscore * fraction * (pident / 100.0)
            
            genus_tid = get_rank_tid(tid, "genus")
            if genus_tid and genus_tid in taxid2name:
                genus_name = taxid2name[genus_tid]
                contig_genus_scores[contig][genus_name] += score
                
                # Update best TaxID tracking for kingdom retrieval
                if contig not in contig_best_genus_tid or score > contig_best_genus_tid[contig][1]:
                    contig_best_genus_tid[contig] = (genus_tid, score)

    # --- 4. Calculate Posterior Probabilities and Export ---
    genus_info = {} # Track kingdom for each genus
    contig_count = defaultdict(int)

    print(f"[*] Exporting posterior probabilities to {posterior_file}...")
    with open(posterior_file, "w") as fout:
        fout.write("contig_id\tbest_genus\tposterior\tkingdom\tall_posteriors\n")
        
        for contig, genus_dict in contig_genus_scores.items():
            total_score = sum(genus_dict.values())
            if total_score == 0: continue
            
            # Calculate probabilities
            genus_probs = {g: s / total_score for g, s in genus_dict.items()}
            best_genus = max(genus_probs, key=genus_probs.get)
            best_post = genus_probs[best_genus]

            # Get Kingdom info
            genus_tid, _ = contig_best_genus_tid.get(contig, (None, 0))
            kingdom_name = "NA"
            if genus_tid:
                sk_tid = get_rank_tid(genus_tid, "superkingdom")
                kingdom_name = taxid2name.get(sk_tid, "NA")

            # Write individual contig results
            sorted_posts = [f"{g}:{p:.3f}" for g, p in sorted(genus_probs.items(), key=lambda x: -x[1])]
            fout.write(f"{contig}\t{best_genus}\t{best_post:.3f}\t{kingdom_name}\t{'; '.join(sorted_posts)}\n")

            # Accumulate for abundance table if posterior threshold (0.8) is met
            if best_post >= 0.8:
                contig_count[best_genus] += 1
                genus_info[best_genus] = kingdom_name

    # --- 5. Export Genus Abundance Table ---
    print(f"[*] Exporting genus abundance to {abundance_file}...")
    with open(abundance_file, "w") as fout:
        fout.write("Genus\tKingdom\tContig_Count\n")
        for genus, count in sorted(contig_count.items(), key=lambda x: x[1], reverse=True):
            fout.write(f"{genus}\t{genus_info.get(genus, 'NA')}\t{count}\n")

if __name__ == "__main__":
    main()
