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

## 🧬 Assembly Comparisons in Figure 3

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
 1. **Assembly length, Number of contigs and Contig N50**
     ```bash
     seqkit stats -a <assembly.fasta>
    ```
 2. **Number of chromosome-level  contigs**
    
    Use the file Chromosome_length.ods and count cells with values ≥99% in columns D, F, and H.
    
 4. **Number of mapped genes**
    ```bash
    minimap2 -cx map-ont <assembly.fasta> ToxoDB-68_TgondiiME49_Genes.fasta > mapping_genes.paf
    awk '{print $1}' mapping_genes.paf | sort | uniq | wc -l
    ```
 5. **Fraction of perfectly mapped genes**
  
    ```bash
    python scripts/Check_the_number_of_genes_with_perfect_matches_in_paf.py mapping_genes.paf
    ```
    
## 📊 Reproducing Figure 3C  

Barplot

    ```bash
    python scripts/comparing_two_assemblies.py <output.pdf>
    ```

## 📊 Reproducing Figure 3D  

This section describes how to explore the structure of rDNA repeats using long reads from the 2015T dataset.
### Step 1: Map 2015T Reads to the Xia et al. Assembly
Use the `non_mit_reads.fastq` file (produced by the filtering pipeline) and map it to the Xia et al. assembly:

   ```bash
   minimap2 -ax map-ont <xia_assembly.fasta> non_mit_reads.fastq \
     | samtools view -bS - \
     | samtools sort -o mapped_sorted.bam
   
   samtools index mapped_sorted.bam
   ```
### Step 2: Extract the Longest Read in the rDNA Region
Target region: JACEHA010000011.1:1,373,482–1,530,609
Run the script to find the longest read mapped to this region:

   ```bash
   python script/find_longes_read_name.py mapped_sorted.bam JACEHA010000011.1 1373482 1530609  
   ```  
This should print the read name: c415285b-0462-4167-b011-16964416454f
Extract the corresponding read from the FASTQ file:
   ```bash
   seqkit grep -p 'c415285b-0462-4167-b011-16964416454f'non_mit_reads.fastq > longest_read.fastq  
   ```
### Step 3: Map Longest Read to rDNA Reference
The reference rDNA sequence is located in Data/rRNA.fasta. Use the following command:
   ```bash
   minimap2 -cx map-ont rRNA.fasta longest_read.fastq > mapping_longest_reads_to_rRNA.paf
   ```
### Step 4: Visualize rDNA Blocks
Generate a visualization of the rDNA block structure on the longest read:
   ```bash
   python script/visualize_rDNA_blocks.py Data/rRNA.fasta mapping_longest_reads_to_rRNA.paf path/to/output/.pdf
   ```
## 📊 Reproducing Figure 4
###  ROP8–ROP2A Region Sequences
The ROP8–ROP2A-related sequences were extracted from the following three assemblies:
- `2015T_assembly.fa`
- `Data/2000B_ROP8_ROP2A_manually_assembled.fa`
- `Data/Xia_et_al_assembly_after_processing_with_companion.fasta.gz`
> **Note:** The last FASTA file corresponds to the Xia et al. assembly, but it was processed with Companion to bring it into the reference orientation.
All extracted sequences are compiled in:  
`Data/ROP8-ROP2A_sequences.fa`
#### Multiple Sequence Alignment (MSA)   
   ```bash
    clustalo -i ROP8-ROP2A_sequences.fa -o clustelo2.fa --outfmt fasta --force --threads 8  
   ```
To visualize the MSA, you can use Jalview. Annotation data is available in `Data/ROP8-ROP2A_sequences.annotation`

#### 📊 Figure 4B

To produce the heatmap, run the following script:

```bash
python script/MSA_analysis_mod.py clustalo2.fa heatmap.pdf
```
To calculate similarity to the reference ROP2A and ROP8 sequences, use:
   ```bash
   python script/similarity_to_reference_ROPs.py  clustelo2.fa similarity.csv
   ```
## 📊 Reproducing Figure 5

The bar plots in Figure 5 are based on the variant analysis pipeline available at:  
🔗 [https://github.com/YomnaGohar/ToxoVar/](https://github.com/YomnaGohar/ToxoVar/)

This pipeline produces a catalog of variants in the following file:  
`{output_dir}/Graph/graph_construction/results/merged_vg_combined_table_placed_ref_for_igv.vcf`

This VCF file was then filtered to:
- Remove variants called as non-reference in 2015T
- Remove loci where the reference allele in 2015T was present at <90% frequency
- Eliminate 6 additional loci that differed before and after graph construction, where manual investigation suggested that the graph-based genotype was likely incorrect or the true underlying genotype remained ambiguous.

The resulting filtered VCF is available at:  
`Data/merged_vg_combined_table_placed_ref_for_igv_no_2015_manual.vcf.gz`

This file can be specified in the `config.yaml` file to continue the remaining steps of the pipeline, which include:
- Stratifying positions by NUMTs, homopolymers, and tandem repeats in `Data/merged_vg_combined_table_placed_ref_for_igv_no_2015_manual.vcf.gz`
- Using Ensembl VEP to predict the functional impact of the variant
  
To stratify the variants before filtering, use the following commands along with the BED files provided in the ToxoVar repository.
```bash
# Extract VCF header
zgrep '^#' merged_vg_combined_table_placed_ref_for_igv.vcf > header.txt
# Variants outside of NUMTs, homopolymers, and tandem repeats
cat header.txt > variants_no_low_complexity.vcf
bedtools intersect -v -a merged_vg_combined_table_placed_ref_for_igv.vcf -b <tandem.bed> <homopolymer.bed> <numts.bed> >> normal.vcf

# Variants in NUMT regions only
cat header.txt > variants_numt.vcf
bedtools intersect -u -a merged_vg_combined_table_placed_ref_for_igv.vcf -b <numts.bed> >> numt.vcf

# Variants in homopolymers (excluding those also in NUMTs)
cat header.txt > temp.vcf
bedtools intersect -u -a merged_vg_combined_table_placed_ref_for_igv.vcf -b <homopolymer.bed> >> temp.vcf
bedtools intersect -v -a temp.vcf -b <numts.bed> > homopolymer.vcf

# Variants in tandem repeats (excluding those in NUMTs and homopolymers)
cat header.txt > temp.vcf
bedtools intersect -u -a merged_vg_combined_table_placed_ref_for_igv.vcf -b <tandem.bed> >> temp.vcf
bedtools intersect -v -a temp.vcf -b <numts.bed> <homopolymer.bed> > tandem.vcf
```



    









