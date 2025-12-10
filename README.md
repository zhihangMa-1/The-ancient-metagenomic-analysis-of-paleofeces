# The ancient metagenomic analysis of paleofeces
This code are analyses that accompanies the manuscript "Human-Animal Interactions in Prehistoric China: Insights from Metagenomic Analysis of Longshan Period Paleofeces"，and allows the reader to replicate the analysis in here. 
## Software and environment dependencies 
All software used for the analysis, including precise version numbers, is listed below. It is **highly recommended** to use a containerized environment (e.g., Docker or Conda) for installation to ensure reproducibility.
### Core bioinformatics tools 
* **FastQC:** v0.11.9
* **BWA:** v0.7.17

## Quality Control and Adapter Trimming
This step uses leeHom (a specialized ancient DNA tool) for adapter trimming and quality filtering. It is followed by filtering sequences shorter than 30 bp using sga.
```bash
# Define working variables
result=/mnt/analysis/mazhihang/Fenbian_analysis/01.data
id=P1 # Example ID for Paleofeces 1
data1=${result}/Paleofeces1_R1.fq.gz
data2=${result}/Paleofeces1_R2.fq.gz
threads=16
min_length=30

# 1.2. Adapter Trimming, QC, and aDNA Damage Flagging (using leeHom)
~/anaconda3/envs/fast/envs/leehom/bin/leeHom \
    -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
    -s AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
    --ancientdna \
    -fq1 $data1 \
    -fq2 $data2 \
    -fqo $id \
    -t $threads

# 1.3. Unzip the Output
gunzip ${id}.fq.gz

# 1.4. Filter Sequences Shorter than 30 bp (using sga)
sga preprocess --dust-threshold=1 -m $min_length ${id}.fq -o ${id}_filtered_trim.fq
```
Generates quality reports (FastQC) and prepares data for plotting the fragment length distribution.
```bash
# Generate basic sequence statistics (using seqkit)
seqkit stats ${id}_filtered_trim.fq > ${id}_filtered_trim.fq.stats
# Generate FastQC report for filtered sequences
fastqc ${id}_filtered_trim.fq -o .
# Extract sequence lengths (using seqkit and awk) and generate fragment length plot (using Rscript)
seqkit fx2tab -l -n -i ${id}_filtered_trim.fq  | awk '{print $2}' > ${i}.length & done
Rscript /home/mazhihang/Script/length_plot.R ${i}.length ${i}_length.pdf ${i} 
```
