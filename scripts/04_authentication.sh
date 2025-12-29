#!/bin/bash
# Description: Targeted mapping to identified animal/plant references and damage profiling
# Usage: bash scripts/04_authentication.sh <input_fastq> <sample_id> <threads> <ref_fasta>
# 1. Parameter Initialization
DATA=$1
ID=$2
THREADS=$3
REF=$4  # Path to the specific animal or plant reference genome
# Define Directory Structure
OUT_DIR="./results/04.Authentication/${ID}"
mkdir -p "${OUT_DIR}"

bwa aln -l 1024 -n 0.01 -t "${THREADS}" "${REF}" "${DATA}" > "${OUT_DIR}/${ID}.sai"
bwa samse -r "@RG\tID:${ID}\tPL:illumina\tSM:${ID}" \
    "${REF}" "${OUT_DIR}/${ID}.sai" "${DATA}" > "${OUT_DIR}/${ID}.sam"
samtools view -Shb -@ "${THREADS}" -q 30 "${OUT_DIR}/${ID}.sam" -o "${OUT_DIR}/${ID}.bam"
samtools sort -@ "${THREADS}" "${OUT_DIR}/${ID}.bam" -o "${OUT_DIR}/${ID}.sort.bam"
samtools index "${OUT_DIR}/${ID}.sort.bam"
dedup -i "${OUT_DIR}/${ID}.sort.bam" -m -o "${OUT_DIR}"
samtools sort "${OUT_DIR}/${ID}.sort_rmdup.bam" -o "${OUT_DIR}/${ID}.rmdup.sort.bam"
samtools index "${OUT_DIR}/${ID}.rmdup.sort.bam"

# Generate mapping depth and coverage statistics via Qualimap
qualimap bamqc -bam "${OUT_DIR}/${ID}.rmdup.sort.bam" -c -outdir "${OUT_DIR}/qualimap_${ID}"
# Calculate deamination patterns at the ends of DNA fragments
mapDamage -i "${OUT_DIR}/${ID}.rmdup.sort.bam" -r "${REF}" -d "${OUT_DIR}/mapdamage_${ID}"

echo "** Authentication for ${ID} against $(basename ${REF}) complete **"
