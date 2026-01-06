#!/usr/bin/env python3
import sys
from collections import defaultdict

# Updated usage to include database paths as arguments
if len(sys.argv) != 7:
    print("Usage: python bayes_hierarchical.py <blast.tsv> <contig_lengths.tsv> <nodes.dmp> <names.dmp> <detailed_out.tsv> <abundance_out.tsv>")
    sys.exit(1)

blast_file = sys.argv[1]
length_file = sys.argv[2]
nodes_file = sys.argv[3]
names_file = sys.argv[4]
detailed_file = sys.argv[5]
abundance_file = sys.argv[6]

TARGET_RANKS = ["species", "genus", "family", "order", "class", "phylum"]
CONFIDENCE_THRESHOLD = 0.8

# === 1. Load Resources ===
contig_len = {}
with open(length_file) as f:
    for line in f:
        cols = line.strip().split()
        if len(cols) >= 2:
            contig_len[cols[0]] = int(cols[1])

parent = {}
rank = {}
with open(nodes_file) as f:
    for line in f:
        cols = line.strip().split("\t|\t")
        tid, p_tid, r = cols[0].strip(), cols[1].strip(), cols[2].strip()
        parent[tid] = p_tid
        rank[tid] = r

taxid2name = {}
with open(names_file) as f:
    for line in f:
        if "scientific name" in line:
            cols = line.strip().split("\t|\t")
            taxid2name[cols[0].strip()] = cols[1].strip()

# Helper function to find node at specific rank
def get_rank_node(tid, target_r):
    curr = tid
    seen = set()
    while curr != "1" and curr in parent and curr not in seen:
        seen.add(curr)
        if rank.get(curr) == target_r:
            return curr
        curr = parent[curr]
    return None

# Helper function to find Superkingdom
def get_superkingdom(tid):
    curr = tid
    seen = set()
    while curr != "1" and curr in parent and curr not in seen:
        seen.add(curr)
        if rank.get(curr) == "superkingdom":
            return taxid2name.get(curr, "Unknown")
        curr = parent[curr]
    return "Unknown"

# === 2. Parse BLAST ===
contig_raw_scores = defaultdict(lambda: defaultdict(float))
with open(blast_file) as f:
    for line in f:
        cols = line.strip().split("\t")
        if len(cols) < 5: continue
        qseqid, staxid = cols[0], cols[1]
        try:
            bitscore, aln_len, pident = float(cols[2]), float(cols[3]), float(cols[4])
        except ValueError:
            continue
            
        if qseqid not in contig_len: continue
        
        frac = min(aln_len / contig_len[qseqid], 1.0)
        score = bitscore * frac * (pident / 100.0)
        contig_raw_scores[qseqid][staxid] += score

# === 3. Classification Logic ===
final_assignments = {} # contig -> (rank, name, posterior, superkingdom)

for contig, taxid_dict in contig_raw_scores.items():
    assigned = False
    best_tid = max(taxid_dict.items(), key=lambda x: x[1])[0]
    skingdom = get_superkingdom(best_tid)

    for r_level in TARGET_RANKS:
        level_scores = defaultdict(float)
        for tid, s in taxid_dict.items():
            r_node = get_rank_node(tid, r_level)
            if r_node:
                level_scores[taxid2name.get(r_node, r_node)] += s
        
        if not level_scores: continue
            
        total_s = sum(level_scores.values())
        best_name, best_s = max(level_scores.items(), key=lambda x: x[1])
        posterior = best_s / total_s
        
        if posterior >= CONFIDENCE_THRESHOLD:
            final_assignments[contig] = (r_level, best_name, posterior, skingdom)
            assigned = True
            break
            
    if not assigned:
        final_assignments[contig] = ("superkingdom", skingdom, 1.0, skingdom)

# === 4. Export Results ===
with open(detailed_file, "w") as f:
    f.write("Contig_ID\tSuperkingdom\tAssigned_Rank\tTaxon_Name\tPosterior\n")
    for c, (lvl, name, post, sk) in final_assignments.items():
        f.write(f"{c}\t{sk}\t{lvl}\t{name}\t{round(post, 4)}\n")

abundance = defaultdict(int)
taxon_to_sk = {}
for c, (lvl, name, post, sk) in final_assignments.items():
    abundance[(lvl, name)] += 1
    taxon_to_sk[(lvl, name)] = sk

with open(abundance_file, "w") as f:
    f.write("Superkingdom\tRank\tTaxon_Name\tContig_Count\n")
    sorted_abundance = sorted(abundance.items(), key=lambda x: x[1], reverse=True)
    for (lvl, name), count in sorted_abundance:
        sk = taxon_to_sk[(lvl, name)]
        f.write(f"{sk}\t{lvl}\t{name}\t{count}\n")
