#!/bin/bash
# scripts/01_preprocessing.sh

ID=$1
DATA1=$2
DATA2=$3
THREADS=$4
MIN_LENGTH=30
OUT_CLEAN="${ID}_clean.fq"

if [ -f "$OUT_CLEAN" ]; then
    echo "[INFO] Pre-processed file '$OUT_CLEAN' already exists."
    echo "[INFO] Verifying file integrity..."
    if [ -s "$OUT_CLEAN" ]; then
        echo "[INFO] Integrity check passed. Skipping trimming and merging steps."
        exit 0
    else
        echo "[WARNING] '$OUT_CLEAN' is empty. Proceeding with raw data processing..."
    fi
fi

echo "[PROCESS] Executing leeHom for $ID..."
leeHom \
    -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
    -s AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
    --ancientdna \
    -fq1 "$DATA1" \
    -fq2 "$DATA2" \
    -fqo "$ID" \
    -t "$THREADS"

gunzip -c "${ID}.fq.gz" > "${ID}.fq"
sga preprocess --dust-threshold=1 -m $MIN_LENGTH "${ID}.fq" -o "$OUT_CLEAN"
rm "${ID}.fq"
