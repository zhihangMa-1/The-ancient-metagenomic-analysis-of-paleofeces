# The ancient metagenomic analysis of paleofeces
This code are analyses that accompanies the manuscript "Human-Animal Interactions in Prehistoric China: Insights from Metagenomic Analysis of Longshan Period Paleofeces"，and allows the reader to replicate the analysis in here. 
## Software and environment dependencies 
All software used for the analysis, including precise version numbers, is listed below. It is **highly recommended** to use a containerized environment (e.g., Docker or Conda) for installation to ensure reproducibility.
### Core bioinformatics tools 
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



## Quality Control and Adapter Trimming
This step uses leeHom (a specialized ancient DNA tool) for adapter trimming and quality filtering. It is followed by filtering sequences shorter than 30 bp using sga.
```bash
# Define working variables
result=/mnt/analysis/mazhihang/Fenbian_analysis/01.data
id=P1 # Example ID for Paleofeces
data1=${data}/Paleofeces1_R1.fq.gz
data2=${data}/Paleofeces1_R2.fq.gz
threads=12
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
## Taxonomic profiling
### 1. Initial Taxonomic Profiling (KrakenUniq)
```bash
# Define working variables (these would be passed when running the script)
# input_fastq: Path to the clean, host-filtered FastQ file (e.g., P1_meta.fq.gz)
# threads: Number of CPU threads to use (e.g., 16)
# sample: Sample ID (e.g., P1)

# Define file paths based on script arguments
input_fastq=$1
threads=$2
results=/mnt/analysis/mazhihang/Fenbian_analysis/09.KrakenUniq
sample=$3
KRAKEN_DB=/mnt/peaks/krakendb # Define database path

# Run KrakenUniq Classification
krakenuniq --db $KRAKEN_DB \
    --fastq-input ${input_fastq} \
    --threads $threads \
    --output $results/${sample}/${sample}_sequences.krakenuniq \
    --report-file $results/${sample}/${sample}_krakenuniq.output \
    --gzip-compressed \
    --only-classified-out
# 
python /home/mazhihang/Script/KrakenUniq/allrank_filter_krakenuniq.py /mnt/store2/mazhihang/FNQZ/${i}/${i}_krakenuniq.output 1000 200
# 添加注释信息
python3 /home/mazhihang/Script/KrakenUniq/get_lineasges_all.py ${i}/${i}_krakenuniq.output.species.filtered $i/krakenuniq.output.species.filtered.with_lineage.tsv
#生成物种丰度表
python /home/mazhihang/Script/KrakenUniq/generate_abundance_matrix_lineasges.py /mnt/store2/mazhihang/FNQZ samplename.txt Species_abundance_matrix_1000_200.csv
```
### 2. Eukaryotic Refinement (Assembly-based)
#### Deduplication & Assembly
```bash
clumpify.sh in=${DATA_DIR}/${ID}_rl.fq out=${RESULT_DIR}/${ID}_dedup.fq dedupe threads=$THREADS qin=33
megahit -r ${RESULT_DIR}/${ID}_dedup.fq --min-contig-len 500 \
    --num-cpu-threads $THREADS --out-dir ${RESULT_DIR}/megahit_${ID} \
    --out-prefix ${ID} --preset meta-large
```
#### Alignment & Damage Authentication
```bash
bowtie2-build ${RESULT_DIR}/megahit_${ID}/${ID}.contigs.fa ${RESULT_DIR}/${ID}_index
bowtie2 -x ${RESULT_DIR}/${ID}_index -U ${DATA_DIR}/${ID}_rl.fq -S ${RESULT_DIR}/${ID}.sam -p ${THREADS} -N 1
samtools view -bS ${RESULT_DIR}/${ID}.sam | samtools sort -o ${RESULT_DIR}/${ID}.sorted.bam
samtools markdup -r -s ${RESULT_DIR}/${ID}.sorted.bam ${RESULT_DIR}/${ID}.rmdup.bam
samtools index ${RESULT_DIR}/${ID}.rmdup.bam
```
#### Run PyDamage

```bash
cd ${RESULT_DIR}
pydamage analyze ${ID}.rmdup.bam -m 500 -p ${THREADS} -pl -f
#### Filter for ancient DNA signatures
awk -F',' 'NR>1 && $2 >= 0.7 && $12 < 0.05 && $15 >= 5 && $17 >= 0.05 {print $1}' pydamage_results/pydamage_results.csv > ancient_contigs.list
seqtk subseq megahit_${ID}/${ID}.contigs.fa ancient_contigs.list > ancient_contigs.fa

```
#### Taxonomic Assignment & Bayesian Filter
```bash
blastn -query ${RESULT_DIR}/ancient_contigs.fa -db ${NT_DB} \
    -evalue 1e-5 -perc_identity 85 -max_target_seqs 100 \
    -outfmt "6 qseqid staxids bitscore length pident evalue stitle" \
    -num_threads ${THREADS} -out ${RESULT_DIR}/nt_${ID}_confirm.tsv
python ${SCRIPT_DIR}/bayes_genus.py ${RESULT_DIR}/nt_${ID}_confirm.tsv ${RESULT_DIR}/ancient_contigs.fa ${RESULT_DIR}/output.tsv
```
