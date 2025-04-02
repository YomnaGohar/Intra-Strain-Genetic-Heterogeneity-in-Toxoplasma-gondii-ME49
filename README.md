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
1. Generate basic read statistics (number of reads, mean length, N50): ``` seqkit stats -a toxo.fastq ```  
2. Count number of long reads >= 10kb:  ``` seqkit seq -m 10000 toxo.fastq | grep '^@' | wc -l ```  
3. Count number of long reads >= 100kb:  ``` seqkit seq -m 100000 toxo.fastq | grep '^@' | wc -l```

