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
### 📁 Assembly Statistics Resources

The `resources/` directory in this repository includes two helpful tables:

- **`Contig_names.ods`**  
  Maps assembled contigs to chromosome names, based on genome annotation output from Companion (`pseudo.pseudochr.agp` file).

- **`Chromosome_length.ods`**  
Provides chromosome lengths across the three assemblies, based on the longest contigs representing each chromosome as listed in Contig_names.ods..
---

## 📊 Reproducing Figure 3A

This section outlines how to compute assembly statistics and gene mapping metrics for the three assemblies.
 1. **Assembly Statistics**
     ```bash
     seqkit stats -a assembly.fasta
    ```









- counting genes and psudogenes
/home/yomna/counting_number_of_genes.py
