#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 02:33:17 2025

@author: yomna
"""
from Bio import AlignIO
import numpy as np
import pandas as pd
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
end_position = 1748
extracted_data, extracted_sequences = extract_msa_region(msa_file, start_position, end_position)
keys_to_remove = [
    "2000B_ROP2A-2B_stitched_read:36681-38986_from_Augustus_only",
    "2000B_ROP2A-2A_stitched_read:28091-30066_from_Augustus_only",
    "2015T_assembly_ROP2A-2B_contig_46:7313278-7315583_from_Augustus_only",
    "Xia_et_al_assembly_ROP2A-2A_JACEHA010000012.1:7310951-7312866_from_Augustus_only",
]

for key in keys_to_remove:
    extracted_sequences.pop(key, None)  
pid_df = calculate_percent_identity(extracted_sequences)
pid_df = pid_df.astype(float)
pid_df = pid_df.combine_first(pid_df.T)
columns_to_keep = ["ROP8_ME_ref", "ROP2A_ME49_ref"]
pid_df = pid_df[columns_to_keep]
pid_df.to_csv( sys.argv[2])




