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
### Initial Taxonomic Profiling (KrakenUniq)
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
### Eukaryotic Refinement (Assembly-based)
```bash
#1.Deduplication & Assembly
clumpify.sh in=${DATA_DIR}/${ID}_rl.fq out=${RESULT_DIR}/${ID}_dedup.fq dedupe threads=$THREADS qin=33
megahit -r ${RESULT_DIR}/${ID}_dedup.fq --min-contig-len 500 \
    --num-cpu-threads $THREADS --out-dir ${RESULT_DIR}/megahit_${ID} \
    --out-prefix ${ID} --preset meta-large
#2.Alignment & Damage Authentication
bowtie2-build ${RESULT_DIR}/megahit_${ID}/${ID}.contigs.fa ${RESULT_DIR}/${ID}_index
bowtie2 -x ${RESULT_DIR}/${ID}_index -U ${DATA_DIR}/${ID}_rl.fq -S ${RESULT_DIR}/${ID}.sam -p ${THREADS} -N 1
samtools view -bS ${RESULT_DIR}/${ID}.sam | samtools sort -o ${RESULT_DIR}/${ID}.sorted.bam
samtools markdup -r -s ${RESULT_DIR}/${ID}.sorted.bam ${RESULT_DIR}/${ID}.rmdup.bam
samtools index ${RESULT_DIR}/${ID}.rmdup.bam
#3.Run PyDamage
cd ${RESULT_DIR}
pydamage analyze ${ID}.rmdup.bam -m 500 -p ${THREADS} -pl -f
#4.Filter for ancient DNA signatures
awk -F',' 'NR>1 && $2 >= 0.7 && $12 < 0.05 && $15 >= 5 && $17 >= 0.05 {print $1}' pydamage_results/pydamage_results.csv > ancient_contigs.list
seqtk subseq megahit_${ID}/${ID}.contigs.fa ancient_contigs.list > ancient_contigs.fa
#5.Taxonomic Assignment & Bayesian Filter
blastn -query ${RESULT_DIR}/ancient_contigs.fa -db ${NT_DB} \
    -evalue 1e-5 -perc_identity 85 -max_target_seqs 100 \
    -outfmt "6 qseqid staxids bitscore length pident evalue stitle" \
    -num_threads ${THREADS} -out ${RESULT_DIR}/nt_${ID}_confirm.tsv
python ${SCRIPT_DIR}/bayes_genus.py ${RESULT_DIR}/nt_${ID}_confirm.tsv ${RESULT_DIR}/ancient_contigs.fa ${RESULT_DIR}/output.tsv

```
## Authentication and Deamination profile

```bash
data=$1
ID=$2
threads=$3
ref=/mnt/analysis/mazhihang/Fenbian_analysis/05.haplotype/1.ref_sequence/NC_002008.4.fa
bwa aln -l 1024 -n 0.01 -t $threads ${ref} ${data} > ${ID}.sai
bwa samse -r "@RG\tID:foo_lane\tPL:illumina\tLB:library\tSM:${ID}" ${ref} ${ID}.sai ${data} > ${ID}.sam
samtools view -Shb -@ $threads ${ID}.sam -q 30 -o ${ID}.bam #MAPQ大于30bp的
samtools sort -@ $threads ${ID}.bam -o ${ID}.sort.bam
samtools index ${ID}.sort.bam
dedup -i ${ID}.sort.bam -m -o .
samtools sort ${ID}.sort_rmdup.bam -o ${ID}.rmdup.sort.bam
samtools index ${ID}.rmdup.sort.bam  && echo "** Dedup Complete **"
qualimap bamqc -bam ${ID}.rmdup.sort.bam -c -outdir ./dedup_{ID} 
mapDamage -i ${ID}.rmdup.sort.bam -r $ref -d ./mapdamage_${ID}

```

## Mitochondrial phylogenetic tree of Canis lupus
### Step 1: Construction of the Genus-level Consensus Reference
To minimize mapping bias, a consensus mitochondrial sequence was constructed from representative Canis species. Species included: Canis lupus, Canis anthus, Canis latrans, Canis simensis, and Canis rufus (Accession numbers in Table S XXX).
```bash
cat *_NCBI_mitogenome_references.fa > cat_NCBI_mitogenome_references.fa
mafft --thread n cat_NCBI_mitogenome_references.fa > Aln_NCBI_mitogenome_references.fa
```
#### Then build consensus sequence for mitogenome reference sequences downloaded from NCBI

1.) Alignment file opened in Geneious, consensus sequences created with 75% Majority rule for family level/each clade

2.) Alignment created from all clade-consensus mitogenome references in Geneious (*_NCBI_mitogenome_references.fa)

3.) Consensus sequence created with 75% Majority rule from all clade-consensus mitogenome references in Geneious (Aln_cons_Canis.fa)

### Step 2: Mapping and Post-processing
```bash
# Index the consensus reference
bwa index Aln_cons_Canis.fa
# Mapping using BWA aln (optimized for short ancient reads) & Filtering: retain MQ >= 25 and remove unmapped reads
ref=/mnt/analysis/mazhihang/Fenbian_analysis/04.mit_tree/Aln_cons_Canis.fa
bwa aln -l 1024 -n 0.001 -t $threads ${ref} ${data_path}/${id}_rl.fq | bwa samse ${ref}  - ${data_path}/${id}_rl.fq  | samtools view -F 4 -q 25 -@ ${threads} -uS - | samtools sort -@ ${threads} -o ${id}.mt.sort.q25.bam
```
### Step 3: Sample Consensus Generation
```bash
angsd -dofasta 2 -docounts 1 -minmapq 25 -minq 25 -uniqueonly 1 -setMinDepthInd 5 -i AS022106_aln2mt.sort.q20.bam -out Cons_Sample.taxa.depth5
gunzip Cons_Sample.taxa.depth5.fa.gz
```
### Step 4: Concatenating and aligning all Canis mitogenome reference sequences and sample consensus
```bash
cat Cons_Sample.taxa.depth5.fa *_NCBI_mitogenome_references.fa > cat_NCBI_mitogenome_references_query.fa
mafft --thread n cat_NCBI_mitogenome_references_query.fa > Aln_NCBI_mitogenome_references_query.fa
```
### Step 5: Running BEAST (Phylogenetic placement mtDNA)
We confirmed the phylogenetic placement of our sequence using a selection of Elephantidae mitochondrial reference sequences, GTR+G, strict clock, a birth-death substitution model, and ran the MCMC chain for 20,000,000 runs, sampling every 20,000 steps. Convergence was assessed using Tracer v1.7.2 and an effective sample size (ESS) > 200.

1.) Aln_NCBI_mitogenome_references_query.fa opened in BEAUti (v1.10.4)

2.) Run with GTR+G, strict clock, a birth-death substitution model, and ran the MCMC chain for 20,000,000 runs, sampling every 20,000 steps

```bash
beast2 -threads n Aln_NCBI_mitogenome_references_query.xml
```
3.) Tracer (v1.7.2) checked for convergence

4.) TreeAnnotator (v1.10.4), 10% Burnin removed, Maximum Clade Credibility Tree created

5.) FigTree (v1.4.4), Tree visualized with posterior probabilities

## Whole genome analysis of Canis lupus familiaris

#### Step 1：Alignment and Initial Processing
The computational processing of the nuclear genome began with the alignment of pre-processed collapsed reads to the Canis lupus familiaris reference genome (UU_Cfam_GSD_1.0) supplemented with the Y chromosome (NC_051844.1).
```bash
ref=/mnt/analysis/mazhihang/my_DB/dog/Canis_lupus_familiaris.fasta
bwa aln -l 1024 -n 0.01 -t $threads ${ref} ${data} > ${ID}.sai
bwa samse -r "@RG\tID:foo_lane\tPL:illumina\tLB:library\tSM:${ID}" ${ref} ${ID}.sai ${data} > ${ID}.sam
samtools view -Shb -@ $threads ${ID}.sam -o ${ID}.bam
samtools sort -@ $threads ${ID}.bam -o ${ID}.sort.bam
samtools index ${ID}.sort.bam
dedup -i ${ID}.sort.bam -m -o .
samtools sort ${ID}.sort_rmdup.bam -o ${ID}.rmdup.sort.bam
samtools index ${ID}.rmdup.sort.bam  
samtools flagstat ${ID}.rmdup.sort.bam > ${ID}.rmdup.sort.bam.flagstat
#qualimap check
qualimap bamqc -bam ${ID}.rmdup.sort.bam -c -outdir ./dedup_{ID} 
mapDamage -i ${ID}.rmdup.sort.bam -r $ref -d ./mapdamage_${ID}
```
#### Step 2：Trim 5bp to remove deamination signals
```bash
bam trimBam ${ID}.mapped.bam ${ID}.trim.bam 5
```
#### Step 3: SNP Calling
```bash
samtools mpileup -R -B -q 30 -Q 30 -l /mnt/rawdata/mazhihang/my_DB/Dog_SNp/Dog10k_Chr.snp.position.filtered  -f /mnt/rawdata/mazhihang/my_DB/Dog_SNp/Canis_lupus_familiars_Chr.fasta ${dir}/1.trimBam/${ID}.trim.bam > pileup_${ID}.txt
pileupCaller --randomHaploid --sampleNames ${ID}  --samplePopName ${ID}  -f /mnt/rawdata/mazhihang/my_DB/Dog_SNp/Dog10k_Chr.snp.filtered -e ${ID} <  pileup_${ID}.txt > pileuplog${ID}.log
```

### Sex Determination
Biological sex for each ancient sample was determined using the Rx ratio method, which calculates the normalized ratio of X-chromosome to autosomal read coverage, a robust approach for low-coverage shotgun sequencing data. The analytical pipeline followed the methodology established by de Flamingh et al. (2020).
```bash
samtools view -q 30 -b input.sam > output.bam
samtools rmdup input.bam output_no_duplicates.bam  
samtools sort input.bam -o sorted_input.bam        
samtools index sorted_input.bam                    
samtools idxstats sorted_input.bam > sorted_input.idxstats
Rscript RX_identifier.R AS24051601 > AS24051601.sex.txt
```

### Genetic relatedness
```bash
#Step 1: Run KINgaroo for data preparation 
KINgaroo \
    -bam ${BAM_DIR} \
    -bed ${BED_FILE} \
    -T ${BAM_LIST} \
    -cnt 0 \
    -c ${THREADS} \
    -o ${KIN_OUT}/kingaroo_results

#Step 2: Run KIN for kinship estimation
KIN -I ${KIN_OUT}/kingaroo_results -O ${KIN_OUT} -c ${THREADS}
```
### Principal Component Analysis (PCA)
```bash
smartpca -p pca.par
```


