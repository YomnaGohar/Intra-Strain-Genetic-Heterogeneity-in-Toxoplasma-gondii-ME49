#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 01:07:51 2025

@author: yomna
"""

import pandas as pd
import sys
import matplotlib.pyplot as plt
gff_file_path = "/home/yomna/hpc_project/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/scaffold.out.gff3"

# Parse GFF file to extract gene descriptions
pseuduo = 0
genes= 0
feature_type_l=[]
with open(gff_file_path, 'rt') as gff:
    for line in gff:
        if line.startswith("#"):
            continue  # Skip comments
        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue  # Skip malformed lines

        chrom, source, feature_type, start, end, score, strand, phase, attributes = cols
        if feature_type == "pseudogene":
           pseuduo +=1
        elif feature_type == "gene":
            genes +=1
        else:
            feature_type_l.append(feature_type)
total=    pseuduo+     genes
print(total) 
print(pseuduo) 
gff_file_path = "/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/Data/t-gondii-me49-new-genome-annotation_results/scaffold.out.gff3"
gff_file_path = "/home/yomna/hpc_project/xie_et_al/xia-et-al_results/scaffold.out.gff3"
# Parse GFF file to extract gene descriptions
pseuduo = 0
genes= 0
with open(gff_file_path, 'rt') as gff:
    for line in gff:
        if line.startswith("#"):
            continue  # Skip comments
        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue  # Skip malformed lines
        chrom, source, feature_type, start, end, score, strand, phase, attributes = cols
        if feature_type == "pseudogene":
           pseuduo +=1
        elif feature_type == "gene":
            genes +=1  
total=    pseuduo+     genes
print(total) 
print(pseuduo)             
