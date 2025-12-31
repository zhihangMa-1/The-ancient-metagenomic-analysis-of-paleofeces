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

#### 1.2 Fragment Length Plot and Quality Assessment (Corresponds to Table S1 and Figure S1)
```bash
#1. Generate comprehensive sequence statistics (Table S1)
seqkit stats "${ID}_clean.fq" > "${ID}_clean_stats.txt"

# 2. Extract fragment lengths
seqkit fx2tab -l -n -i "${ID}_clean.fq" | awk '{print $2}' > "${ID}.lengths"

# 3. Plot length distribution (Figure S1)
# Script: scripts/visualize_length.R
# This plot provides evidence of aDNA degradation by showing short fragment enrichment.
Rscript scripts/visualize_length.R "${ID}.lengths" "${ID}_length_plot.pdf" "${ID}"
```
### 2. Metagenomic Taxonomic Profiling
This section describes the taxonomic identification of ancient metagenomic reads and the generation of the abundance matrix (Corresponds to **Table S3**, **Figure 2**, **Figure 7**).
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
* Key Output: Abundance_matrix.csv (The primary data source for **Figure 2** and **Figure 7**).


### 2.2 Assembly-based Metagenomic Validation and Authentication

To further validate the taxonomic assignments and confirm the ancient origin of identified taxa, we performed de novo assembly followed by damage-based authentication. This step ensures that the sequences used for taxonomic refinement exhibit typical ancient DNA post-mortem deamination patterns.

#### Workflow Summary
* Purpose: Reconstruct longer contigs from reads, authenticate their ancient status via deamination patterns, and refine taxonomic assignments.
* Core Script: scripts/03_blastn.sh (An integrated wrapper for assembly and authentication).
* Key Evidence: This process supports the authenticity of the taxa presented in **Figure 6** and provides the data for the taxonomic refinements.

#### Step 1: Metagenomic Assembly and Read Mapping
We utilize BBTools clumpify for sequence deduplication to improve assembly contiguity, followed by MEGAHIT for de novo assembly. Reads are then mapped back to the resulting contigs to evaluate damage.
* Usage：
```bash
# Arguments: <sample_id> <input_fastq> <threads> <nt_db_path>
bash scripts/03_blastn.sh sample_01 data/clean/sample01_clean.fq 16 /path/to/ncbi_nt_db
```
#### Step 2: Ancient DNA Damage Authentication
* Filtering Criteria for "Ancient" Status:

    * Accuracy: >= 0.7 (Ensures high confidence in the damage model).

    * Q-value: < 0.05 (Statistical significance of the deamination signal).

    * Deamination Signal: Must be observed at fragment termini.

* Primary Output: ancient_contigs.fa (A FASTA file containing only authenticated ancient sequences).

#### Step 3:Taxonomic Assignment and Bayesian Refinement
Authenticated contigs are queried against the NCBI nt database using blastn. A custom Bayesian model is then applied to resolve potential taxonomic ambiguities between closely related genera.
* Script: scripts/bayes_genus.py
* Command:
```bash
python scripts/bayes_genus.py \
    results/sample_01/nt_sample_01_confirm.tsv \
    results/sample_01/ancient_contigs.fa \
    results/sample_01/sample_01_final_taxonomic_refinement.tsv
```

### 2.3. Targeted Authentication and Deamination Profile
This section describes the targeted mapping of reads to specific animal or plant reference genomes to evaluate DNA damage patterns. This step is essential for authenticating the ancient origin of the taxa identified in the metagenomic screening.
#### Workflow Overview
* Purpose: Verify the authenticity of specific taxa by analyzing post-mortem DNA damage (C-to-T deamination).
* Key Script: scripts/04_authentication.sh
* Direct Evidence:

    * **Figure S2,S3,S4 and S5**: DNA damage plots (generated by mapDamage).
    * **Table S4**: Mapping statistics including depth and coverage (generated by Qualimap).

We use BWA aln with parameters optimized for ancient DNA (stringent gap opening and seeding) to map reads against specific references (e.g., Dog, Pig, or dietary plants).After removing duplicates, we quantify mapping quality and ancient DNA signatures.
* Usage:
  
```bash
# Arguments: <input_fastq> <sample_id> <threads> <ref_fasta>
bash scripts/04_authentication.sh data/clean/sample01.fq sample01 16 refs/dog_genome.fa
```
* Key Outputs:
    * Mapping Stats: Located in ${OUT_DIR}/qualimap_${ID} (Source for **Table S4**).
    * Damage Plots: Located in ${OUT_DIR}/mapdamage_${ID}/Fragmisincorporation_plot.pdf (Source for **Figure S2,S3,S4**).

## 3. Mitochondrial phylogenetic tree of Canis lupus

This section details the reconstruction of the mitochondrial genome (mtDNA) and the subsequent Bayesian phylogenetic inference used to determine the maternal lineage of the identified Canis lupus.
### Step 1: Genus-level Consensus Reference Construction
To avoid "reference bias" toward a single modern dog or wolf genome, we constructed a consensus reference from multiple Canis species (Accession numbers in Table S1).
```bash
cat *_NCBI_mitogenome_references.fa > cat_NCBI_mitogenome_references.fa
mafft --thread n cat_NCBI_mitogenome_references.fa > Aln_NCBI_mitogenome_references.fa
```
### Step 2: Then build consensus sequence for mitogenome reference sequences downloaded from NCBI

* Alignment file opened in Geneious, consensus sequences created with 75% Majority rule for family level/each clade

* Alignment created from all clade-consensus mitogenome references in Geneious

* Consensus sequence created with 75% Majority rule from all clade-consensus mitogenome references in Geneious 

### Step 3: Mapping and Consensus Generation

* Usage:
```bash
# Arguments: <sample_id> <input_fastq> <threads> <consensus_ref>
bash scripts/05_mtDNA_consensus.sh sample_01 data/clean.fq 16 refs/canis_consensus.fa
```
* Key Parameters & Tools:

    * BWA aln: -l 1024 -n 0.001 (Optimized for ancient DNA to capture short, divergent fragments).
    * ANGSD: -setMinDepthInd 5 (Requires at least 5-fold coverage per site to ensure consensus accuracy).

* Primary Output: ${ID}_mt_consensus.fa (The reconstructed mitochondrial sequence).

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
* Final Result: The MCC tree corresponds to **Figure 3 (Phylogenetic tree of Canis mtDNA)**.

## 4. Paleogenomic Analysis of Canis lupus familiaris

This section describes the genomic pipeline for analyzing dog paleogenomes, including mapping, sex determination, SNP calling, and population structure analysis.

#### 4.1 Alignment and Post-processing
Reads were aligned to the Canis lupus familiaris reference genome (UU_Cfam_GSD_1.0) supplemented with the Y chromosome (NC_051844.1). We utilized a specialized ancient DNA pipeline to mitigate the impact of post-mortem deamination.

* Workflow:
    * Alignment: bwa aln -l 1024 -n 0.01 (stringent mapping for ancient DNA).
    * Deduplication: Performed via dedup to remove PCR/optical duplicates.
    * Deamination Trimming: The first 5 bp of both ends were trimmed using bam trimBam to remove damaged bases before downstream SNP calling.

* Usage:
```bash
# Arguments: <input_fastq> <sample_id> <threads> <ref_fasta>
bash scripts/06_genome_processing.sh sample_01 data/clean.fq 16 dog_ref.fa
```
#### 4.2 Sex Determination (Corresponds to Table S5)

Biological sex was determined using the Rx ratio method (de Flamingh et al., 2020), which calculates the ratio of X-chromosome coverage to the average autosomal coverage.

```bash
# Step 1: Generate coverage statistics
samtools idxstats sample_01.final.bam > sample_01.idxstats
# Step 2: Run R script for Rx ratio calculation
Rscript scripts/RX_identifier.R "sample_01" "sample_01.idxstats"
```
* Results: Outputs the biological sex (Male/Female) and confidence intervals (Corresponds to **Table S5**).
  
#### 4.3: SNP Calling

To minimize bias in low-coverage ancient genomes, we performed pseudo-haploid SNP calling against the Dog10K dataset.

* Software: samtools mpileup and pileupCaller.

* Key Parameters:

    * samtools mpileup -q 30 -Q 30: High-quality filters for mapping and base quality.

    * pileupCaller --randomHaploid: Randomly selects one allele per site to account for low coverage.

* Input Files: SNP positions are defined in configs/Dog10k_SNPs.bed.

* Output: EIGENSTRAT format files (.geno, .snp, .ind) for population genetics.

### 4.4 Kinship and Population Structure

To identify biological relationships and genetic affinities, we utilized KIN and smartPCA.

#### Kinship Estimation (Corresponds to **Table S5** )

Biological relationships between samples were estimated using KIN.

* Usage:  
```bash
# Step 1: Data preparation
KINgaroo -bam "${OUT_DIR}" -bed "$SNP_POS" -T configs/bam.list -cnt 0 -c "${THREADS}" -o "${OUT_DIR}/kinship_prep"
# Step 2: Estimate kinship coefficients
KIN -I "${OUT_DIR}/kinship_prep" -O "${OUT_DIR}/kinship_results" -c "${THREADS}"
```

#### Principal Component Analysis (PCA) (Corresponds to **Figure 4** )
Population structure was visualized by projecting ancient samples onto modern dog reference panel.
```bash
# Run smartPCA with lsqproject: YES
smartpca -p configs/pca.par > "${OUT_DIR}/pca.log"
```
* Key Result: This analysis generates the data for **Figure 4 (PCA plot)**.

