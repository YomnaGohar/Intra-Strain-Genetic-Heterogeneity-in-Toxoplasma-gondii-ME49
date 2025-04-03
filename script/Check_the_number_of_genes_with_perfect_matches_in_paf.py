#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 17:13:16 2024

@author: yomna
"""

import pandas as pd
import sys
paf_file_1 = sys.argv[1]
def fraction_perfectly_mapped_genes(paf_file):
    # Define columns to parse from PAF
    columns_to_use = [0, 1, 2, 3, 5, 6, 7, 8, 9,10]
    paf = pd.read_csv(
        paf_file,
        sep="\t",
        header=None,
        usecols=columns_to_use,
        names=["query_name", "query_length", "query_start", "query_end",
               "target_name", "target_length", "target_start", "target_end",
               "match_length","match_length_gap"]
    )

    # Calculate query coverage and alignment length
    paf["query_coverage"] = (paf["query_end"] - paf["query_start"]) / paf["query_length"]

    # Identify perfectly mapped genes
    paf["perfect_mapping"] = (
       # (paf["query_coverage"] == 1.0) &  # Full query coverage (no clipping)
        (paf["match_length"] == paf["query_length"]) &  # No mismatches or insertions
        (paf["match_length_gap"] == paf["query_length"])  # No deletions
    )

    # Count perfectly mapped genes
    perfectly_mapped_genes = paf[paf["perfect_mapping"]]["query_name"].nunique()

    # Total number of unique genes
    total_genes = paf["query_name"].nunique()

    # Fraction of perfectly mapped genes
    fraction = perfectly_mapped_genes / total_genes if total_genes > 0 else 0

    return fraction, perfectly_mapped_genes, total_genes

# Count genes with high query coverage in each assembly
genes_with_high_coverage_my_assembly = fraction_perfectly_mapped_genes(paf_file_1)
#genes_with_high_coverage_xie_assembly = fraction_perfectly_mapped_genes(paf_file_xia)

# Print the results
print(f"Number of genes with  perfect match: {genes_with_high_coverage_my_assembly}")
