#!/bin/bash
# Description: Read merging and adapter trimming for aDNA
# Usage: bash scripts/01_preprocessing.sh <ID> <R1_fastq> <R2_fastq> <threads>

ID=$1
DATA1=$2
DATA2=$3
THREADS=$4
MIN_LENGTH=30

# 1. Adapter Trimming and Read Merging via leeHom
# Using universal Illumina adapter sequences
leeHom \
    -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
    -s AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
    --ancientdna \
    -fq1 "$DATA1" \
    -fq2 "$DATA2" \
    -fqo "$ID" \
    -t "$THREADS"

# 2. Filtering sequences shorter than 30bp and removing low-complexity reads (Dust threshold = 1)
gunzip -c "${ID}.fq.gz" > "${ID}.fq"
sga preprocess --dust-threshold=1 -m $MIN_LENGTH "${ID}.fq" -o "${ID}_clean.fq"
