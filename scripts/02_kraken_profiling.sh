#!/bin/bash
# Description: Taxonomic Profiling via KrakenUniq
# Usage: bash scripts/02_kraken_profiling.sh <input_fastq> <threads> <sample_id> <db_path>

INPUT_FASTQ=$1
THREADS=$2
SAMPLE=$3
KRAKEN_DB=$4
RESULTS_DIR=$5

mkdir -p "$RESULTS_DIR/$SAMPLE"

# 1. Run KrakenUniq Classification
krakenuniq --db "$KRAKEN_DB" \
    --fastq-input "$INPUT_FASTQ" \
    --threads "$THREADS" \
    --output "$RESULTS_DIR/$SAMPLE/${SAMPLE}_sequences.krakenuniq" \
    --report-file "$RESULTS_DIR/$SAMPLE/${SAMPLE}_krakenuniq.output" \
    --gzip-compressed \
    --only-classified-out
# 2. Dual-threshold Filtering (1000 unique k-mers, 200 reads)
python scripts/allrank_filter_krakenuniq.py \
    "$RESULTS_DIR/$SAMPLE/${SAMPLE}_krakenuniq.output" 1000 200
# 3. Lineage Annotation
python3 scripts/get_lineasges_all.py \
    "$RESULTS_DIR/$SAMPLE/${SAMPLE}_krakenuniq.output.species.filtered" \
    "$RESULTS_DIR/$SAMPLE/krakenuniq.output.species.filtered.with_lineage.tsv"
# 4. Abundance Matrix Generation
python scripts/generate_abundance_matrix_lineasges.py \
    "$RESULTS_DIR" samplename.txt Abundance_matrix.csv
