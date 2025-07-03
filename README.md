#  PRO-seq, precision nuclear run-on and sequencing #

For everyone who has a new face to HaoLab, we make this pipeline/protocol to standardize upstream analysis pipeline for PRO-seq data.
Any questions or trouble u have met during this procudure please conntact me, AKA. Liu Qian. Here I post my domestic E-mail, to reach me in time: [18229048546@163.com](mailto:18229048546@163.com)
## Part I Introduction ##
### i. Workflow ###
Here stands an throughout workflow of PRO-seq data analysis.

![](https://github.com/jieniqianqian/PRO-seq/tree/main/Figure/PRO-seq.drawio.png)

As illustrated in the figure,
(i) yellow circles represent the steps where commands need to be entered;
(ii) pink dashed rectangular boxes represent the output results after processing at each step.

We will proceed by structuring our workflow according to (ii).
### ii. File Structure ###
Here stands an throughout file structure of PRO-seq data analysis.

- *You can decide what your structure looks like, which makes it more efficient to work.*
```
PRO-seq 
├─ 1.raw_data 
│    └─ QC 
├─ 2.clean_data
├─ 3.chr
├─ 4.bam
|	 ├─ uniq_bam
|	 ├─ fail_qc_dups_bam
│    └─ noMT_bam
├─ 5.bw
│    ├─ plus_bw
│    └─ minus_bw
├─ 6.bed
|    ├─ pause_index
|    └─ mRNA_contamination
└─ 7.stat
       ├─ pause_index
       ├─ mRNA_contamination
       └─ TSS_enrichment
```
### iii. mamba Environment ###
You can configure a mamba environment named 'PRO-seq' using the following code, which includes the essential software for PRO-Seq analysis.

```bash
micromamba create -n PRO-seq falco fastp seqtk sortmerna bowtie2 samtools deeptools
```
## Part II Generation of Data for Analysis: FASTQ2BAM ##
In this section, you will convert the raw FASTQ files into BAM files, which can be used for subsequent analysis.

### i. Raw Data Quality Check(QC) ###

You can perform quality check on the raw data to assess the sequencing quality.

```bash
falco -o PRO-seq/1.raw_data/QC/ PRO-seq/1.raw_data/raw_data_1.fq.gz PRO-seq/1.raw_data/raw_data_2.fq.gz
```
### ii. Trimming After QC ###
Based on the quality control results, you can perform appropriate trimming on the raw data.

```bash
fastp \
--overrepresentation_analysis \
--thread 16 \
--in1 PRO-seq/1.raw_data/raw_data_1.fq.gz  --adapter_sequence AGATCGGAAGAG \
--trim_poly_g \
--low_complexity_filter \
--length_required 20 \
--average_qual 20 \
-h PRO-seq/2.clean_data/1.html \
-j PRO-seq/2.clean_data/1.json \
-o PRO-seq/2.clean_data/1.noadap.fq

seqtk trimfq \
-b 1 PRO-seq/2.clean_data/1.noadap.fq \ 
> PRO-seq/2.clean_data/1.trimmed.fq

seqtk seq \
-L 20 \
-r PRO-seq/2.clean_data/1.trimmed.fq \
> PRO-seq/2.clean_data/1.temp.fq

sed -e 's|\\:[^:]*\\([[:space:]].*\\)|\\1 |g' PRO-seq/2.clean_data/1.temp.fq  > PRO-seq/2.clean_data/1.clean.fq
```
### iii. Mapping Clean Data to rDNA first to eliminating rRNA ###
Align the cleaned data obtained in the previous step to the reference rDNA. After eliminating the abundant rRNA,the resulting data will be prepared for further analysis.

```bash
sortmerna --ref rfam-5.8s-database-id98.fasta \
          --ref rfam-5s-database-id98.fasta \
          --ref silva-euk-18s-id95.fasta \
          --ref silva-euk-28s-id98.fasta \
          --reads PRO-seq/2.clean_data/1.clean.fq \
          --aligned PRO-seq/2.clean_data/1.align \
          --other PRO-seq/2.clean_data/1.non_rRNA \
          --workdir PRO-seq/2.clean_data/ \
          --fastx \
          --log -v \
          -a 20 \
          --out2 T
```
### iV. Mapping Clean non rRNA Data to Genome ###

Align the cleaned data obtained in the previous step to the reference genome of the corresponding species. After filtering for high-quality sequences, the resulting data will be prepared for further analysis.
#### 1. Mapping, Filtering ####
Align the quality-controlled FASTQ files to the reference genome, gets rid of multimapped reads, and perform sorting.

```bash
bowtie2 -p 16 \
--very-sensitive \
-q --end-to-end \
--un-conc PRO-seq/2.clean_data/1.hgRemovedTemp.fastq \
--rg-id R1 \
-x reference_species_index \
-U PRO-seq/2.clean_data/1.non_rRNA.fq | \
awk '$2 != 4 {print}' | \
sed '/XS:/d' | \
samtools view -bS - -@ 1 | \
samtools sort - -@ 1 -o PRO-seq/4.bam/uniq_bam/uniq.sorted.bam
```

#### 2. BAM Indexing ####
Index the sorted files.

```
samtools index PRO-seq/4.bam/uniq_bam/uniq.sorted.bam
```

#### 3. Filter for high-quality aligned sequences  ####

Split genome mapping result bamfile into two: high-quality aligned reads (keepers) and unmapped reads.

```bash
samtools view -q 30 -b -@ 16 -U PRO-seq/4.bam/fail_qc_dups_bam/fail_qc_dups.bam PRO-seq/4.bam/uniq_bam/uniq.sorted.bam > PRO-seq/4.bam/uniq_bam/high_quality_sort.bam
```
#### 4. BAM Indexing ####

Index the high quality sorted files.

```bash
samtools index PRO-seq/4.bam/uniq_bam/high_quality_sort.bam
```
### V. Remove mitochondrial reads  ###

Remove reads that are mapped to mitochondria.

#### 1. Acquire non-mitochondrial chromosomes ####  

Exactly match non-mitochondrial chromosomes from high quality sorted files, where chrM and chrMT are mitochondria in UCSC, M and MT are mitochondria in Ensembl, and rCRSd is the revised Cambridge mitochondria.

```bash
samtools idxstats PRO-seq/4.bam/high_quality_sort.bam | \
cut -f 1-2 | awk '{print $1, 0, $2}' | \
grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > PRO-seq/3.chr/chr_sizes.bed
```
#### 2. Obtain reads on non-mitochondrial chromosomes ####  

Remove mitochondrial reads based on the non-mitochondrial chromosomes.

```bash
samtools view -L PRO-seq/3.chr/chr_sizes.bed -b -@ 16 PRO-seq/4.bam/uniq_bam/high_quality_sort.bam > PRO-seq/4.bam/noMT_bam/noMT.bam
```
#### 3. BAM Indexing ####

Index the non-mitochondrial high quality sorted files.

```bash
samtools index PRO-seq/4.bam/noMT_bam/noMT.bam
```
### VI. TSS enrichment  ###

For Assessment of  run-on efficiency, aggregates sequencing reads at 2000 bases upstream and downstream of a reference set of TSSs to plot and calculate a TSS enrichment score. The normalized TSS enrichment score is the ratio of the average coverage in 100-bp windows, with the numerator centered at the TSS peak summit and the denominator in the background at the edge of the 2000-bp window.

#### 1. Extract chromosome names and their lengths ####

```bash
samtools view -H PRO-seq/4.bam/noMT_bam/noMT.bam | \
grep 'SN:' | \
awk -F':' '{print $2,$3}' | \
awk -F' ' -v OFS='\t' '{print $1,$3}' > PRO-seq/3.chr/chr_order.txt
cut -f 1 PRO-seq/3.chr/chr_order.txt > PRO-seq/3.chr/chr_keep.txt
```
#### 2. Plus TSS enrichment ####

Using Python scripts from PEPPRO(https://github.com/databio/peppro/tree/master/tools).

```bash
python pyTssEnrichment.py \
-a PRO-seq/4.bam/noMT_bam/noMT.bam -b reference_plus_TSS -p ends \
-c 16 \
-z -v -s 6 -o PRO-seq/7.stat/TSS_enrichment/plus_TssEnrichment.txt
```
#### 3. Minus  TSS enrichment ####

Using Python scripts from PEPPRO(https://github.com/databio/peppro/tree/master/tools).

```bash
python pyTssEnrichment.py \
-a PRO-seq/4.bam/noMT_bam/noMT.bam -b reference_minus_TSS -p ends \
-c 16 \
-z -v -s 6 -o PRO-seq/7.stat/TSS_enrichment/minus_TssEnrichment.txt
```
#### 4. Plot  plus TSS Enrichment  ####

Using R scripts from PEPPRO(https://github.com/databio/peppro/tree/master/tools).

```bash
Rscript PEPPRO.R tss -i PRO-seq/7.stat/TSS_enrichment/plus_TssEnrichment.txt
```
#### 5. Plot  minus TSS Enrichment  ####

Using R scripts from PEPPRO(https://github.com/databio/peppro/tree/master/tools).

```bash
Rscript PEPPRO.R tss -i PRO-seq/7.stat/TSS_enrichment/minus_TssEnrichment.txt
```
### VII. Pause index  ###

Define the pause index as the ratio of the density of reads in the *pausing* region versus the density in the corresponding gene body,then plots the frequency distribution of the pause index across genes. A greater pause index indicates a more efficient run-on, as a higher value indicates that paused polymerases efficiently incorporate the modified NTPs. An efficient run-on process has a median pause index greater than 10.

#### 1. Remove missing chr from TSS and gene body annotations ####

Sort reference TSSes and reference gene bodies.

```bash
grep -wf PRO-seq/3.chr/chr_keep.txt reference_pi_tss | bedtools sort -i stdin -faidx PRO-seq/3.chr/chr_order.txt > PRO-seq/3.chr/chr_keep.txt > PRO-seq/6.bed/pause_index/pi_tss.bed
grep -wf PRO-seq/3.chr/chr_keep.txt reference_pi_body | bedtools sort -i stdin -faidx PRO-seq/3.chr/chr_order.txt > PRO-seq/3.chr/chr_keep.txt > PRO-seq/6.bed/pause_index/pi_body.bed
```

#### 2. Determine coverage of highest scoring TSS ####

Sort TSS annotations by read count.

```bash
bedtools coverage -sorted -counts -s -a PRO-seq/6.bed/pause_index/pi_tss.bed -b PRO-seq/4.bam/noMT_bam/noMT.bam -g PRO-seq/3.chr/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > PRO-seq/6.bed/pause_index/TSS_density.bed
```
#### 3. Determine coverage of highest scoring gene body ####

Sort gene body annotations by read count.

```bash
bedtools coverage -sorted -counts -s -a PRO-seq/6.bed/pause_index/pi_body.bed -b PRO-seq/4.bam/noMT_bam/noMT.bam -g PRO-seq/3.chr/chr_order.txt | awk '$7>0' | sort -k4 > PRO-seq/6.bed/pause_index/gene_body_density.bed
```
#### 4. Deduplicate Gene  ####
```bash
cat PRO-seq/6.bed/pause_index/TSS_density.bed PRO-seq/6.bed/pause_index/gene_body_density.bed | awk '{print $4}' > PRO-seq/7.stat/pause_index/PI_all_genes.txt
awk -F';' '$1 in first{print first[$1] $0; first[$1]=""; next} {first[$1]=$0 ORS}' PRO-seq/7.stat/pause_index/PI_all_genes.txt | awk '!NF {print;next}; !($0 in a) {a[$0];print}' > PRO-seq/7.stat/pause_index/PI_shared_genes.txt
```

#### 5. Find genes in TSS annotation  ####
```bash
awk -F '\t' 'NR==FNR {id[$1]; next} $4 in id' PRO-seq/7.stat/pause_index/PI_shared_genes.txt PRO-seq/6.bed/pause_index/TSS_density.bed > PRO-seq/6.bed/pause_index/TSS_shared_density.bed
```
#### 6. Find genes in gene body annotation  ####
```bash
awk -F '\t' 'NR==FNR {id[$1]; next} $4 in id' PRO-seq/7.stat/pause_index/PI_shared_genes.txt PRO-seq/6.bed/pause_index/gene_body_density.bed > PRO-seq/6.bed/pause_index/gene_body_shared_density.bed
```

#### 7. Calculate the background pause index #### 

Get the background value of TSS and gene body of the same gene.

```bash
awk 'BEGIN{FS=OFS="\t"} FNR>0 && FNR==NR{a[$4]=$4 OFS $0; next} FNR>0{print $0,a[$4]?a[$4]:"\t"}' PRO-seq/6.bed/pause_index/TSS_shared_density.bed PRO-seq/6.bed/pause_index/gene_body_shared_density.bed | \
awk -v OFS='\t' '{ if ($6 == "+"){print $9, $10, $3, $4,sqrt((($15+$7)/sqrt(($3-$10)^2))^2),($15/sqrt(($11-$10)^2))/($7/sqrt(($3-$2)^2)), $6} else {print $9, $10, $3, $12,sqrt((($15+$7)/sqrt(($10-$2)^2))^2),($15/sqrt(($11-$10)^2))/($7/sqrt(($3-$2)^2)), $6}}' | \
env LC_COLLATE=C sort -k1,1 -k2,2n > PRO-seq/6.bed/pause_index/temp.bed
```

#### 8. Filter background value #### 

Only values above the median of the pause index background value are valid.

```bash
median_expr=$(awk '{print $5}' PRO-seq/6.bed/pause_index/temp.bed | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}')
awk -v OFS='\t' -v median_expr=$median_expr '{ if ($5 > median_expr) {print $1, $2, $3, $4, $6, $7}}' PRO-seq/6.bed/pause_index/temp.bed > PRO-seq/6.bed/pause_index/pause_index.bed
```
#### 9. Plot  pause index   ####

Using R scripts from PEPPRO(https://github.com/databio/peppro/tree/master/tools).

```bash
Rscript PEPPRO.R pi --annotate -i PRO-seq/6.bed/pause_index/pause_index.bed
```

### VIII. mRNA contamination  ###

For Assessment of  nascent RNA purity,calculates the exon to intron read density ratio. And ChRO-seq libraries have an increased promoter emphasis and higher mRNA contamination indicated by an increase in reads in promoters and exons at the cost of reads in introns and promoter flanking regions .

#### 1. Remove missing chr from exon and intron annotations ####

Sort reference exons and reference introns.

```bash
grep -wf PRO-seq/3.chr/chr_keep.txt reference_exon | bedtools sort -i stdin -faidx PRO-seq/3.chr/chr_order.txt > PRO-seq/3.chr/chr_keep.txt > PRO-seq/6.bed/mRNA_contamination/exons_sort.bed
grep -wf PRO-seq/3.chr/chr_keep.txt reference_intron | bedtools sort -i stdin -faidx PRO-seq/3.chr/chr_order.txt > PRO-seq/3.chr/chr_keep.txt >  PRO-seq/6.bed/mRNA_contamination/introns_sort.bed
```

#### 2. Determine coverage of highest scoring exons ####

Sort exon annotations by read count.

```bash
bedtools coverage -sorted -counts -s -a PRO-seq/6.bed/mRNA_contamination/exons_sort.bed -b PRO-seq/4.bam/noMT_bam/noMT.bam -g PRO-seq/3.chr/chr_order.txt > PRO-seq/6.bed/mRNA_contamination/exons_coverage.bed
```
#### 3. Determine coverage of highest scoring introns ####

Sort intron annotations by read count.

```bash
bedtools coverage -sorted -counts -s -a PRO-seq/6.bed/mRNA_contamination/introns_sort.bed -b PRO-seq/4.bam/noMT_bam/noMT.bam -g PRO-seq/3.chr/chr_order.txt > PRO-seq/6.bed/mRNA_contamination/introns_coverage.bed
```
#### 4. Determine exonic and intronic RPKM for individual genes  ####

Need Total Reads divided by 1M.

```bash
ar=$(samtools view -F 4 -c PRO-seq/4.bam/noMT_bam/noMT.bam)
scaling_factor=$(echo "scale=6; $ar/1000000" | bc)
awk -v OFS='\t' -v scaling_factor=$scaling_factor '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/scaling_factor)/geneSizeKB[a], strand[a] }}' PRO-seq/6.bed/mRNA_contamination/exons_coverage.bed | \
awk '$5>0' | sort -k4 > PRO-seq/6.bed/mRNA_contamination/exons_rpkm.bed
awk -v OFS='\t' -v scaling_factor=$scaling_factor '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/scaling_factor)/geneSizeKB[a], strand[a] }}' PRO-seq/6.bed/mRNA_contamination/introns_coverage.bed | \
awk '$5>0' | sort -k4 > PRO-seq/6.bed/mRNA_contamination/introns_rpkm.bed
```
#### 5. Deduplicate Gene  ####
```bash
cat PRO-seq/6.bed/mRNA_contamination/exons_rpkm.bed PRO-seq/6.bed/mRNA_contamination/introns_rpkm.bed | awk '{print $4}' > PRO-seq/7.stat/mRNA_contamination/all_genes.txt
awk -F';' '$1 in first{print first[$1] $0; first[$1]=""; next} {first[$1]=$0 ORS}' PRO-seq/7.stat/mRNA_contamination/all_genes.txt | awk '!NF {print;next}; !($0 in a) {a[$0];print}' > PRO-seq/7.stat/mRNA_contamination/shared_genes.txt
```

#### 6. Find genes in exonic RPKM #### 

```bash
awk -F '\t' 'NR==FNR {id[$1]; next} $4 in id' PRO-seq/7.stat/mRNA_contamination/shared_genes.txt PRO-seq/6.bed/mRNA_contamination/exons_rpkm.bed > PRO-seq/6.bed/mRNA_contamination/shared_exons.bed
```

#### 7. Find genes in intronic RPKM #### 

```bash
awk -F '\t' 'NR==FNR {id[$1]; next} $4 in id' PRO-seq/7.stat/mRNA_contamination/shared_genes.txt PRO-seq/6.bed/mRNA_contamination/introns_rpkm.bed > PRO-seq/6.bed/mRNA_contamination/shared_introns.bed
```
#### 8. Calculate the mRNA contamination   ####

```bash
awk 'BEGIN{FS=OFS="\t"} FNR>0 && FNR==NR{a[$4]=$4 OFS $0; next} FNR>0{print $0,a[$4]?a[$4]:"\t"}' PRO-seq/6.bed/mRNA_contamination/shared_exons.bed PRO-seq/6.bed/mRNA_contamination/shared_introns.bed | \
awk 'BEGIN { OFS="\t" } {print $8, $9, $10, $11, ($12/$5), $13}' | sort -k1,1 -k2,2n > PRO-seq/6.bed/mRNA_contamination/exon_intron_ratios.bed
```

#### 9. Plot  mRNA contamination   ####

Using R scripts from PEPPRO(https://github.com/databio/peppro/tree/master/tools).

```bash
Rscript PEPPRO.R mrna -i PRO-seq/6.bed/mRNA_contamination/exon_intron_ratios.bed --annotate
```
### IX. Peak calling  ###

Identify regions of read enrichment based on the alignment results.

```bash
bamCoverage -b PRO-seq/4.bam/noMT_bam/noMT.bam -o bw/plus_bw/plus.bw --samFlagExclude 16 -bs 50 --normalizeUsing CPM -p 16
bamCoverage -b PRO-seq/4.bam/noMT_bam/noMT.bam -o bw/minus_bw/minus.bw --samFlagInclude 16 -bs 50 --normalizeUsing CPM -p 16
```

### OUTPUT ###

After completing the above steps, you will obtain 4 figures similar to the ones shown in 'example_report.ppt', slide (1-4),which can be used to evaluate the quality of the PRO-seq library.

#### OUTPUT 1 ####

standard:main peak>10(Log10 main peak > 1).

#### OUTPUT 2 ####

standard:main peak 1-1.8(0 <Log10 main peak < 0.26).

#### OUTPUT 3 ####

standard:main peak 25-50 bp downstream of the TSS.

