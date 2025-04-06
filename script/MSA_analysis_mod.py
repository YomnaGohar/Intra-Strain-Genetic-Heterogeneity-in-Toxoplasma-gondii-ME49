#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 02:33:17 2025

@author: yomna
"""
from Bio import AlignIO
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
def extract_msa_region(msa_file, start, end):
    alignment = AlignIO.read(msa_file, "fasta")  # Load MSA
    extracted_data = {}
    extracted_sequences = {}
    print(f"Extracting bases from positions {start} to {end} and computing lengths & gap positions...\n")
    for record in alignment:
        sequence_id = record.id
        sequence = str(record.seq)
        extracted_seq = sequence[start-1:end]  # Convert to 0-based indexing
        ungapped_length = len(extracted_seq.replace("-", ""))  # Remove gaps and calculate length
        gap_positions = [i + start for i, base in enumerate(extracted_seq) if base == "-"]
        extracted_data[sequence_id] = {
            "ungapped_length": ungapped_length,
            "gap_positions": gap_positions
        }
        extracted_sequences[sequence_id] = extracted_seq
    return extracted_data, extracted_sequences

def calculate_percent_identity(sequences):
    seq_ids = list(sequences.keys())
    num_seqs = len(seq_ids)    
    pid_matrix = np.full((num_seqs, num_seqs), np.nan)
    for i in range(num_seqs):
        for j in range(i, num_seqs):  
            seq1 = sequences[seq_ids[i]]
            seq2 = sequences[seq_ids[j]]
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != "-")
            total_valid_positions = sum(1 for a, b in zip(seq1, seq2) if not (a == "-" and b == "-"))            
            pid = (matches / total_valid_positions) * 100 if total_valid_positions > 0 else 0
            pid_matrix[i, j] = pid  # Fill matrix at both (i, j) and (j, i)
    pid_df = pd.DataFrame(pid_matrix, index=seq_ids, columns=seq_ids)
    return pid_df



msa_file= sys.argv[1] 
start_position = 409
end_position = 2212
extracted_data, extracted_sequences = extract_msa_region(msa_file, start_position, end_position)
keys_to_remove = [
    "ROP2A_ME49_ref",
    "2000B_ROP2A-2B_stitched_read:36681-38986_from_Augustus_only",
    "2000B_ROP2A-2A_stitched_read:28091-30066_from_Augustus_only",
    "2015T_assembly_ROP2A-2B_contig_46:7313278-7315583_from_Augustus_only",
    "Xia_et_al_assembly_ROP2A-2A_JACEHA010000012.1:7310951-7312866_from_Augustus_only",
    "ROP8_ME_ref"
]

for key in keys_to_remove:
    extracted_sequences.pop(key, None)  # safely removes the key if it exists

pid_df = calculate_percent_identity(extracted_sequences)
pid_df = pid_df.astype(float)
pid_df = pid_df.combine_first(pid_df.T)

np.fill_diagonal(pid_df.values, 100)

# Set the size of the figure
plt.figure(figsize=(15, 12))  # Increase figure size

# Create a clustered heatmap with "rocket" palette and larger annotation text
g = sns.clustermap(
    pid_df, 
    method="average",  # Clustering method
    metric="euclidean",  # Distance metric
    cmap="rocket",  # Use Rocket colormap
    linewidths=0.5, 
    figsize=(15, 12),  # Increase figure size
    annot=True,  # Display numbers inside each cell
    fmt=".1f",  # Format numbers with one decimal place
    annot_kws={"size": 10},  # Increase font size of numbers
    cbar_kws={'label': 'Percent Identity'}
)

plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=12)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=12)
output_image = sys.arg[2] 
plt.savefig(output_image, dpi=300, bbox_inches="tight")

