# Human-Animal Interactions in Prehistoric China: Metagenomic Analysis of Longshan Period Paleofeces
This repository contains the computational pipeline and scripts to accompany the manuscript "Human-Animal Interactions in Prehistoric China: Insights from Metagenomic Analysis of Longshan Period Paleofeces" by Ma et al., submitted to the Journal of Archaeological Science. This compendium provides all necessary code, configuration files, and documentation required to reproduce the taxonomic profiling, authentication, and paleogenomic analyses presented in the study. The repository is permanently archived and accessible via DOI: XXX.

## Repository Structure
This repository centralizes all scripts and configurations required to replicate the study:

* scripts/: Core analysis pipeline (Bash, Python, and R scripts).
* configs/: Configuration templates, SNP lists, and example metadata.
* LICENSE: MIT License.
* README.md: Main documentation.

## Software and Environment 
All software used for the analysis, including precise version numbers, is listed below. We have provided a Docker image [Link] that encapsulates all necessary software dependencies and precise tool versions.

### Core Tools & Versions
* Processing: FastQC (v0.11.9), leeHom, sga, seqkit.
* Taxonomy: KrakenUniq (v0.5.8), MEGAHIT (v1.2.9), blastn.
* Ancient DNA: pydamage (v0.70), mapDamage (v2.2.1).
* Genomics: BWA (v0.7.17), samtools (v1.13), angsd (v0.935), pileupCaller.
* Population Genetics: BEAST2 (v2.6.6), smartPCA (EIG v7.2.1), KIN (v1.0).


## Detailed Analysis Workflow

### 1. Data Pre-processing
#### 1.1 Adapter trimming and Filtering

Process raw reads including adapter trimming and length filtering.
```bash
# Usage: bash scripts/01_preprocessing.sh <sample_id> <fq1> <fq2> <threads>
bash scripts/01_preprocessing.sh sample_01 sample_01_R1.fq.gz sample_01_R2.fq.gz 8
```

#### 1.2 Fragment Length Plot and Quality Assessment
```bash
#1. Generate comprehensive sequence statistics (Corresponds to Table S1)
seqkit stats "${ID}_clean.fq" > "${ID}_clean_stats.txt"

# 2. Extract fragment lengths
seqkit fx2tab -l -n -i "${ID}_clean.fq" | awk '{print $2}' > "${ID}.lengths"

# 3. Plot length distribution (Corresponds to Figure S1)
# Script: scripts/visualize_length.R
# This plot provides evidence of aDNA degradation by showing short fragment enrichment.
Rscript scripts/visualize_length.R "${ID}.lengths" "${ID}_length_plot.pdf" "${ID}"
```
### 2. Metagenomic Taxonomic Profiling
This section describes the taxonomic identification of ancient metagenomic reads and the generation of the abundance matrix (Corresponds to Table S3, Figure 2, Figure 7).
#### 2.1 K-mer Profiling and dual-threshold Filtering

##### Step 1: Run Taxonomic Classification

* Script: scripts/02_kraken_profiling.sh
* Usage:
```bash
# Arguments: <input_fastq> <threads> <sample_id> <db_path> <results_dir>
bash scripts/02_kraken_profiling.sh \
    data/clean/sample01_clean.fq 16 sample01 \
    /path/to/kraken_db ./results/taxa
```
* Inputs: Pre-processed clean FASTQ files.
* Outputs: ${SAMPLE}_krakenuniq.output: Comprehensive k-mer report. ${SAMPLE}_sequences.krakenuniq: Per-read classification results.

##### Step 2: Dual-threshold Filtering & Lineage Annotation
We provide a Python-based pipeline to filter the KrakenUniq output and assign full taxonomic lineages.
* Scripts: scripts/allrank_filter_krakenuniq.py and scripts/get_lineasges_all.py
* Usage:
```bash
# 1. Filter results (Thresholds: 1000 unique k-mers, 200 reads)
python scripts/allrank_filter_krakenuniq.py \
    results/taxa/sample01/sample01_krakenuniq.output 1000 200

# 2. Annotate lineages for the filtered species list
python scripts/get_lineasges_all.py \
    results/taxa/sample01/sample01_krakenuniq.output.species.filtered \
    results/taxa/sample01/taxa_with_lineage.tsv
```

##### Step 3: Abundance Matrix Generation (Corresponds to Table S3)
This script aggregates filtered results from all samples into a single matrix for downstream statistical analysis and visualization in R.
* Script: scripts/generate_abundance_matrix_lineasges.py
* Usage:
```bash
# Arguments: <results_directory> <sample_list_file> <output_csv_name>
python scripts/generate_abundance_matrix_lineasges.py \
    ./results/taxa scripts/samplename.txt Abundance_matrix.csv
```
* Key Output: Abundance_matrix.csv (The primary data source for Figure 2 and Figure 7).


### 2.2 Assembly-based Metagenomic Validation and Authentication

To further validate the taxonomic assignments and confirm the ancient origin of identified taxa, we performed de novo assembly followed by damage-based authentication. This step ensures that the sequences used for taxonomic refinement exhibit typical ancient DNA post-mortem deamination patterns.

#### Workflow Summary
* Purpose: Reconstruct longer contigs from reads, authenticate their ancient status via deamination patterns, and refine taxonomic assignments.
* Core Script: scripts/03_blastn.sh (An integrated wrapper for assembly and authentication).
* Key Evidence: This process supports the authenticity of the taxa presented in Figure 6 and provides the data for the taxonomic refinements.

#### Step 1: Metagenomic Assembly and Read Mapping
We utilize BBTools clumpify for sequence deduplication to improve assembly contiguity, followed by MEGAHIT for de novo assembly. Reads are then mapped back to the resulting contigs to evaluate damage.
* Usage：
```bash
# Arguments: <sample_id> <input_fastq> <threads> <nt_db_path>
bash scripts/03_blastn.sh sample_01 data/clean/sample01_clean.fq 16 /path/to/ncbi_nt_db
```
#### Step 2: Ancient DNA Damage Authentication
We employ pyDamage to analyze the C-to-T deamination frequency at the ends of DNA fragments for each contig.
    * Filtering Criteria for "Ancient" Status:
        * Accuracy: $\ge$ 0.7
        * Q-value: $<$ 0.05Deamination signal: Observed at fragment termini.
    * Primary Output: ancient_contigs.fa (A FASTA file containing only authenticated ancient sequences).


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


