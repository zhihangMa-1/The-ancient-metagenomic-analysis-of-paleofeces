#!/bin/bash
# Description: WGS mapping, deduplication, and deamination trimming
# Usage: bash scripts/06_genome_processing.sh <input_fastq> <sample_id> <threads> <ref_fasta>

DATA=$1
ID=$2
THREADS=$3
REF=$4  # UU_Cfam_GSD_1.0 + Y chromosome

OUT_DIR="./results/06.WGS/${ID}"
mkdir -p "${OUT_DIR}"

bwa aln -l 1024 -n 0.01 -t "${THREADS}" "${REF}" "${DATA}" > "${OUT_DIR}/${ID}.sai"
bwa samse -r "@RG\tID:${ID}\tPL:illumina\tLB:lib1\tSM:${ID}" \
    "${REF}" "${OUT_DIR}/${ID}.sai" "${DATA}" > "${OUT_DIR}/${ID}.sam"
samtools view -Shb -@ "${THREADS}" -q 30 "${OUT_DIR}/${ID}.sam" | \
samtools sort -@ "${THREADS}" -o "${OUT_DIR}/${ID}.sort.bam"
dedup -i "${OUT_DIR}/${ID}.sort.bam" -m -o "${OUT_DIR}"
samtools sort "${OUT_DIR}/${ID}.sort_rmdup.bam" -o "${OUT_DIR}/${ID}.rmdup.sort.bam"
samtools index "${OUT_DIR}/${ID}.rmdup.sort.bam"
qualimap bamqc -bam "${OUT_DIR}/${ID}.rmdup.sort.bam" -c -outdir "${OUT_DIR}/qualimap_${ID}"
mapDamage -i "${OUT_DIR}/${ID}.rmdup.sort.bam" -r "${REF}" -d "${OUT_DIR}/mapdamage_${ID}"
bam trimBam "${OUT_DIR}/${ID}.rmdup.sort.bam" "${OUT_DIR}/${ID}.trim5.bam" 5
samtools sort "${OUT_DIR}/${ID}.trim5.bam" -o "${OUT_DIR}/${ID}.final.bam"
samtools index "${OUT_DIR}/${ID}.final.bam"
