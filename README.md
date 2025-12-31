# Human-Animal Interactions in Prehistoric China: Metagenomic Analysis of Longshan Period Paleofeces
This repository contains the computational pipeline and scripts to accompany the manuscript "Human-Animal Interactions in Prehistoric China: Insights from Metagenomic Analysis of Longshan Period Paleofeces" by Ma et al., submitted to the Journal of Archaeological Science. This compendium provides all necessary code, configuration files, and documentation required to reproduce the taxonomic profiling, authentication, and paleogenomic analyses presented in the study. The repository is permanently archived and accessible via DOI: XXX.

## 1.Repository Structure
This repository centralizes all scripts and configurations required to replicate the study:

* scripts/: Core analysis pipeline (Bash, Python, and R scripts).
* configs/: Configuration templates, SNP lists, and example metadata.
* LICENSE: MIT License.
* README.md: Main documentation.

## 2.Software and Environment 
All software used for the analysis, including precise version numbers, is listed below. We have provided a Docker image [Link] that encapsulates all necessary software dependencies and precise tool versions.

### Core Tools & Versions
* **FastQC:** v0.11.9
* **BWA:** v0.7.17
* **leeHom**
* **sga**
* **seqkit**
* **krakenuniq**
* **megahit**
* **blastn**
* **pydamage**
* **mapDamage**
* **bowtie**
* **samtools**
* **dedup**
* **bwa**
* **mafft**
* **angsd**
* **PileupCaller**
* **KIN**
* **smartPCA**

## 1. Data pre-processing

### 1.1 Read Processing and Adapter Trimming
```bash
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
```

### 1.2 Fragment Length Plot and Quality Assessment
```bash
# Generate comprehensive sequence statistics 
seqkit stats "${ID}_clean.fq" > "${ID}_clean_stats.txt"

# Run FastQC for base quality distribution
fastqc "${ID}_clean.fq" -o ./qc_reports/

# Plot length distribution (R script provided in scripts/visualize_length.R)
seqkit fx2tab -l -n -i "${ID}_clean.fq" | awk '{print $2}' > "${ID}.lengths"
Rscript scripts/visualize_length.R "${ID}.lengths" "${ID}_length_plot.pdf" "${ID}"
```
## 2. Metagenomic Taxonomic Profiling

### 2.1 K-mer Profiling via KrakenUniq

The initial screening used KrakenUniq to count unique k-mers associated with each taxon. To account for potential false positives in ancient samples, we implemented a dual-threshold filter: a species is only considered present if it possesses at least 1,000 unique k-mers and is supported by at least 200 reads.

```bash
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
```
### 2.2 Assembly-based Metagenomic Validation and Authentication

```bash
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
```

## 3. Authentication and Deamination profile

```bash
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
```

## 4. Mitochondrial phylogenetic tree of Canis lupus

### Step 1: Genus-level Consensus Reference Construction
To avoid "reference bias" toward a single modern dog or wolf genome, we constructed a consensus reference from multiple Canis species (Accession numbers in Table S1).
```bash
cat *_NCBI_mitogenome_references.fa > cat_NCBI_mitogenome_references.fa
mafft --thread n cat_NCBI_mitogenome_references.fa > Aln_NCBI_mitogenome_references.fa
```
#### Then build consensus sequence for mitogenome reference sequences downloaded from NCBI

1.) Alignment file opened in Geneious, consensus sequences created with 75% Majority rule for family level/each clade

2.) Alignment created from all clade-consensus mitogenome references in Geneious

3.) Consensus sequence created with 75% Majority rule from all clade-consensus mitogenome references in Geneious 

### Step 2: Mapping and Consensus Generation

```bash
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
```

### Step 4: Running BEAST (Phylogenetic placement mtDNA)
The reconstructed mitochondrial genome was aligned with a global dataset of Canis references to perform Bayesian phylogenetic inference.

```bash
# 1. Concatenate sample consensus with reference dataset
cat "${OUT_DIR}/${ID}_mt_consensus.fa" references/Global_Canis_Dataset.fa > "${OUT_DIR}/${ID}_query_aln.fa"

# 2. Multiple Sequence Alignment
mafft --thread "${THREADS}" "${OUT_DIR}/${ID}_query_aln.fa" > "${OUT_DIR}/${ID}_final_alignment.fa"

# 3. Bayesian Inference (BEAST2)
# Parameters: GTR+G substitution model, Strict Clock, Birth-Death Prior
# MCMC: 20,000,000 iterations, sampling every 20,000 steps
beast2 -threads "${THREADS}" configs/mtDNA_analysis.xml

# 4. Post-MCMC Diagnostics
(1) Convergence: ESS (Effective Sample Size) values were verified using Tracer v1.7.2 (all ESS > 200).
(2) Tree Summarization: A Maximum Clade Credibility (MCC) tree was generated using TreeAnnotator (10% burn-in).
(3) Visualization: Final trees were visualized in FigTree v1.4.4.
```
## 5. Paleogenomic Analysis of Canis lupus familiaris

#### 5.1 Alignment and Post-processing
Reads were aligned to the Canis lupus familiaris reference genome (UU_Cfam_GSD_1.0) supplemented with the Y chromosome (NC_051844.1). We utilized a specialized ancient DNA pipeline to mitigate the impact of post-mortem deamination.

```bash
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
```
#### 5.2 Sex Determination
Biological sex was determined using the Rx ratio method (de Flamingh et al., 2020), which calculates the ratio of X-chromosome coverage to the average autosomal coverage.
```bash
# Calculate alignment statistics
samtools idxstats "${OUT_DIR}/${ID}.final.bam" > "${OUT_DIR}/${ID}.idxstats"

# Identify sex using R script
Rscript scripts/RX_identifier.R "${ID}" "${OUT_DIR}/${ID}.idxstats" > "${OUT_DIR}/${ID}.sex.txt"
```
#### 5.3: SNP Calling
We performed SNP calling against the Dog10K dataset using pileupCaller.

```bash
SNP_POS="./configs/Dog10k_SNPs.bed"
SNP_REF="./configs/Dog10k_SNPs.snp"
REF_FASTA="./references/Canis_lupus_familiaris.fasta"

# Step 1: Generate mpileup at specific SNP locations
# -q 30 -Q 30: Minimum mapping and base quality thresholds
samtools mpileup -R -B -q 30 -Q 30 \
    -l "$SNP_POS" \
    -f "$REF_FASTA" \
    "${OUT_DIR}/${ID}.final.bam" > "${OUT_DIR}/${ID}.mpileup"

# Step 2: Call pseudo-haploid genotypes using pileupCaller
# --randomHaploid: Randomly selects one base per site
# -e: Output in EIGENSTRAT format (geno, snp, ind)
pileupCaller --randomHaploid \
    --sampleNames "${ID}" \
    --samplePopName "${ID}" \
    -f "$SNP_REF" \
    -e "${OUT_DIR}/${ID}" \
    < "${OUT_DIR}/${ID}.mpileup" \
    > "${OUT_DIR}/${ID}.pileup.log"
```

### 5.4 Kinship and Population Structure
To identify biological relationships and genetic affinities, we utilized KIN and smartPCA.

#### Kinship Estimation

```bash
# Step 1: Data preparation
KINgaroo -bam "${OUT_DIR}" -bed "$SNP_POS" -T configs/bam.list -cnt 0 -c "${THREADS}" -o "${OUT_DIR}/kinship_prep"
# Step 2: Estimate kinship coefficients
KIN -I "${OUT_DIR}/kinship_prep" -O "${OUT_DIR}/kinship_results" -c "${THREADS}"
```
#### Principal Component Analysis (PCA)
```bash
# Run smartPCA with lsqproject: YES
smartpca -p configs/pca.par > "${OUT_DIR}/pca.log"
```


