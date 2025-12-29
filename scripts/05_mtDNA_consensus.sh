#!/bin/bash
# Description: mtDNA mapping and sample consensus generation
# Usage: bash scripts/05_mtDNA_consensus.sh <input_fastq> <sample_id> <threads> <consensus_ref>

ID=$1
DATA=$2
THREADS=$3
REF=$4
OUT_DIR="./results/05.mtDNA/${ID}"

mkdir -p "${OUT_DIR}"

# 1. Indexing the reference
bwa index "${REF}"

# 2. Competitive mapping optimized for aDNA
# Using -n 0.001 and -l 1024 for high sensitivity
bwa aln -l 1024 -n 0.001 -t "${THREADS}" "${REF}" "${DATA}" | \
bwa samse "${REF}" - "${DATA}" | \
samtools view -F 4 -q 25 -@ "${THREADS}" -uS - | \
samtools sort -@ "${THREADS}" -o "${OUT_DIR}/${ID}.mt.sort.q25.bam"

# 3. Consensus generation using ANGSD
# We used -setMinDepthInd 5 to ensure base-calling accuracy
angsd -dofasta 2 -docounts 1 -minmapq 25 -minq 25 -uniqueonly 1 \
      -setMinDepthInd 5 -i "${OUT_DIR}/${ID}.mt.sort.q25.bam" \
      -out "${OUT_DIR}/${ID}_mt_consensus"

# Decompress the generated FASTA
gunzip "${OUT_DIR}/${ID}_mt_consensus.fa.gz"
