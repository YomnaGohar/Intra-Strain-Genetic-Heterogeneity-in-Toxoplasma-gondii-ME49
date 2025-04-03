# Intra-Strain Genetic Heterogeneity in *Toxoplasma gondii* ME49

This repository provides a guide to reproducing the core analysis in our study of genomic divergence across independently maintained *T. gondii* ME49 isolates.

---

## 📦 Data Access

Download raw Nanopore sequencing data from **NCBI BioProject**:

- **BioProject:** [PRJNA1241696](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1241696)

### Biosample IDs:
- **2015T:**
- **2020T:** `SAMN47561284`, `SAMN47561285`
- **2000B:** `SAMN47757524` to `SAMN47757542`

---

## 🧹 Data Filtering

To filter raw reads (remove host and mitochondrial DNA), use the pipeline available at:

🔗 [YomnaGohar/T.-gondii_filter_assemble](https://github.com/YomnaGohar/T.-gondii_filter_assemble)

This pipeline outputs for each sample:
- `toxo.fastq` – filtered reads likely from *T. gondii*
- `non_mit_reads.fastq` – *T. gondii* reads excluding mitochondrial DNA

---

## 📊 Reproducing Figure 2B

This section outlines how to compute basic sequencing statistics for each sample.
1. Generate basic read statistics (T. gondii sequencing data, number of reads, mean length, N50): ``` seqkit stats -a toxo.fastq ```  
2. Count number of long reads >= 10kb:  ``` seqkit seq -m 10000 toxo.fastq | grep '^@' | wc -l ```  
3. Count number of long reads >= 100kb:  ``` seqkit seq -m 100000 toxo.fastq | grep '^@' | wc -l```
4. Coverage: T. gondii sequencing data / 65000000

   ---

## 🧬 Assembly Comparisons

This section describes how to prepare the data needed to compare gene mapping across different *Toxoplasma gondii* assemblies for Figure 3A.

---

### 📥 Download Required Assemblies

Download the following genome assemblies from GenBank or RefSeq:

| Assembly | Description | Accession |
|----------|-------------|-----------|
| **2015T** | Custom assembly used in this study | GenBank accession: *[to be added]* |
| **Xia et al.** | ME49 isolate assembly | GenBank: `JACEHA000000000.1` |
| **TGA4** | Reference strain | RefSeq: `GCF_000006565.2` |

---

### 🧬 Prepare Gene FASTA for Mapping

Since ToxoDB does not provide a direct gene FASTA file, you must construct it using the genome and annotation files.

#### 1. Download genome and annotation from [ToxoDB](https://toxodb.org/toxo/app/downloads):
- `ToxoDB-68_TgondiiME49_Genome.fasta`
- `ToxoDB-68_TgondiiME49.gff`

Save both files in your working directory.

#### 2. Generate a BED file with gene coordinates:
```bash
awk -F'\t' '!/Parent=/' ToxoDB-68_TgondiiME49.gff | grep -v '^#' \
| awk '{match($9, /ID=([^;]+)/, arr); print $1, $4-1, $5, arr[1]}' OFS='\t' \
> ToxoDB-68_TgondiiME49.bed
```
#### 3. Extract gene sequences from the genome FASTA:
```bash
bedtools getfasta \
  -fi ToxoDB-68_TgondiiME49_Genome.fasta \
  -bed ToxoDB-68_TgondiiME49.bed \
  -fo ToxoDB-68_TgondiiME49_Genes.fasta \
  -name
```
---

## 📊 Reproducing Figure 3A

This section outlines how to compute asembly statistics for the three assemblies.
1. Assembly length, Number of contigs and Contig N50: ``` seqkit stats -a <assembly.fasta> ```  
2. data directory in this repos contains tables Contig_names.ods and Chromosome_length.ods. the first one was obtained from genome annotain results from Companion pseudo.pseudochr.agp file that maps the assembled contigs to Chromsome name. the number of Chromsoome level contigs was obtained from the largest contig making the chromsoe and divided by the expected length of the chromsome which is the average of the chromsome length in the thress assembles. percentage >=99% is consdered chromsome-level assembly contig
3.  - figure 3A, number of mapped genes



the fasta file is in  /home/yomna/Desktop/PhD_Yomna_Gohar/papers/Toxo_variants_paper_prefinal/Data/ToxoDB-68_TgondiiME49_Genes_2.fasta and contains 8778


conda activate snakemake-tutorial
genes were mapped to xie et al: minimap2 -cx map-ont /home/yomna/hpc/fasta/ToxoME49_GenBank_JACEHA000000000.1.fasta /home/yomna/Desktop/PhD_Yomna_Gohar/papers/Toxo_variants_paper_prefinal/Data/ToxoDB-68_TgondiiME49_Genes_2.fasta > mapping_genes.paf
results are in: /home/yomna/hpc_project/xie_et_al/mapped_genes_to_the_assembly
counting number of mapped genes: awk '{print $1}' mapping_genes.paf | sort | uniq | wc -l

TGA4:
conda activate snakemake-tutorial
minimap2 -cx map-ont GCF_000006565.2_TGA4_genomic.fasta /home/yomna/Desktop/PhD_Yomna_Gohar/papers/Toxo_variants_paper_prefinal/Data/ToxoDB-68_TgondiiME49_Genes_2.fasta >  mapping_genes_to_TGA4_using_new_fasta.paf
results are in: /home/yomna/Desktop/PhD_Yomna_Gohar/Toxo_referance
counting: awk '{print $1}' mapping_genes_to_TGA4_using_new_fasta.paf | sort | uniq | wc -l
python /home/yomna/Desktop/PhD_Yomna_Gohar/papers/Toxo_variants_paper/additional_scripts/Check_the_number_of_genes_with_perfect_matches_in_paf.py  mapping_genes_to_TGA4_using_new_fasta.paf

2015T assembly:
conda activate snakemake-tutorial
minimap2 -cx map-ont filtered_contigs_without_contig_51.fasta /home/yomna/Desktop/PhD_Yomna_Gohar/papers/Toxo_variants_paper_prefinal/Data/ToxoDB-68_TgondiiME49_Genes_2.fasta >  mapping_genes.paf
results are in: /home/yomna/hpc_project/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/40-polishing
counting: awk '{print $1}' mapping_genes.paf | sort | uniq | wc -l
python /home/yomna/Desktop/PhD_Yomna_Gohar/papers/Toxo_variants_paper/additional_scripts/Check_the_number_of_genes_with_perfect_matches_in_paf.py mapping_genes.paf

- counting genes and psudogenes
/home/yomna/counting_number_of_genes.py
