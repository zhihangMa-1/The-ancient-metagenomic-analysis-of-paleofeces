#!/bin/bash
# Description: contig assembly, authentication, and BLAST validation
# Usage: bash scripts/03_blastn.sh <sample_id> <input_fastq> <threads> <nt_db_path>

ID=$1
INPUT_FQ=$2
THREADS=$3
NT_DB=$4

# Define Directory Structure (Generalized Paths)
RESULT_DIR="./${ID}"
mkdir -p "${RESULT_DIR}"

# 1. Assembly
# BBTools clumpify removes optical/PCR duplicates to improve assembly contiguity
clumpify.sh \
    in="${INPUT_FQ}" \
    out="${RESULT_DIR}/${ID}_dedup.fq" \
    dedupe threads="${THREADS}" qin=33
# MEGAHIT assembly with parameters optimized for metagenomes
megahit \
    -r "${RESULT_DIR}/${ID}_dedup.fq" \
    --min-contig-len 500 \
    --num-cpu-threads "${THREADS}" \
    --out-dir "${RESULT_DIR}/megahit_${ID}" \
    --out-prefix "${ID}" \
    --preset meta-large

# 2. Damage Authentication
bowtie2-build "${RESULT_DIR}/megahit_${ID}/${ID}.contigs.fa" "${RESULT_DIR}/${ID}_index"
bowtie2 -x "${RESULT_DIR}/${ID}_index" \
    -U "${RESULT_DIR}/${ID}_dedup.fq" \
    -S "${RESULT_DIR}/${ID}.sam" \
    -p "${THREADS}" -N 1
samtools view -bS "${RESULT_DIR}/${ID}.sam" | samtools sort -o "${RESULT_DIR}/${ID}.sorted.bam"
samtools markdup -r -s "${RESULT_DIR}/${ID}.sorted.bam" "${RESULT_DIR}/${ID}_rmdup.bam"
samtools index "${RESULT_DIR}/${ID}_rmdup.bam"

#3. Ancient DNA Signature Quantification
pydamage analyze "${RESULT_DIR}/${ID}_rmdup.bam" -m 500 -p "${THREADS}" -pl -f

#4. Extract Authenticated Ancient Contigs
# Filtering criteria: accuracy >= 0.7, q-value < 0.05, and deamination signal present
awk -F',' 'NR>1 && $2 >= 0.7 && $12 < 0.05 && $15 >= 5 && $17 >= 0.05 {print $1}' \
    "${RESULT_DIR}/pydamage_results/pydamage_results.csv" > "${RESULT_DIR}/ancient_contigs.list"
# Extract sequences of authenticated contigs
seqtk subseq "${RESULT_DIR}/megahit_${ID}/${ID}.contigs.fa" \
    "${RESULT_DIR}/ancient_contigs.list" > "${RESULT_DIR}/ancient_contigs.fa"

#5. Taxonomic Assignment & Bayesian Filtering
# Query authenticated contigs against the NCBI nt database
blastn -query "${RESULT_DIR}/ancient_contigs.fa" \
    -db "${NT_DB}" \
    -evalue 1e-5 -perc_identity 85 -max_target_seqs 100 \
    -outfmt "6 qseqid staxids bitscore length pident evalue stitle" \
    -num_threads "${THREADS}" \
    -out "${RESULT_DIR}/nt_${ID}_confirm.tsv"

# Resolve taxonomic ambiguities using a custom Bayesian model
python scripts/bayes_genus.py \
    "${RESULT_DIR}/nt_${ID}_confirm.tsv" \
    "${RESULT_DIR}/ancient_contigs.fa" \
    "${RESULT_DIR}/${ID}_final_taxonomic_refinement.tsv"

echo "** Assembly-based Refinement for ${ID} Complete **"
