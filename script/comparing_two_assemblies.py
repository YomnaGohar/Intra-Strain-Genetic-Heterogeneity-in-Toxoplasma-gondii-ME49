#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 13:47:46 2024

@author: yomna
"""
import matplotlib.pyplot as plt
import sys
xie_et_al={
    "JACEHA010000001.1":[1997081],
    "JACEHA010000002.1":[1982890],
    "JACEHA010000003.1":[2497824],
    "JACEHA010000004.1,JACEHA010000032.1":[2591259,7546],
    "JACEHA010000005.1":[2830561],
    "JACEHA010000006.1":[3596493],
    "JACEHA010000007.1":[3781463],
    "JACEHA010000008.1,JACEHA010000009.1":[3410322,1395699],
    "JACEHA010000010.1":[12088238],
    "JACEHA010000011.1":[6603527],
    "JACEHA010000012.1":[7667134],
    "JACEHA010000013.1,JACEHA010000033.1":[6675137,3707],
    "JACEHA010000014.1":[7423918]
    }
isolate_2015T={
    "contig_14":[2000958],
    "contig_37":[1987479],
    "contig_13,contig_22":[2450531,53715],
    "contig_24":[2596606],
    "contig_2":[2881782],
    "contig_4,contig_5":[3514252,45802],
    "contig_20":[3788807],
    "contig_42,contig_40":[4595037,76924],
    "contig_43":[12101396],
    "contig_35,contig_30":[1390229,5219375],
    "contig_46":[7662873],
    "contig_39":[6685531],
    "contig_11,contig_9":[7421671,18981]
    }

# Function to create a combined barplot with X-axis in Mb and no legends, mapping based on order with accumulated bar heights
def create_combined_barplot(data1, data2):
    fig, ax = plt.subplots(figsize=(10, 6))

    # List of contigs from both assemblies (order matters, so use the same order)
    contigs1 = list(data1.keys())
    contigs2 = list(data2.keys())
    
    # Ensure both lists have the same length and are aligned
    if len(contigs1) != len(contigs2):
        raise ValueError("Both assemblies must have the same number of contigs")

    # Loop through the contigs and plot data for both assemblies
    for i, (contig1, contig2) in enumerate(zip(contigs1, contigs2)):
        sizes1 = data1[contig1]
        sizes2 = data2[contig2]
        
        # Start plotting Xie et al. Assembly
        if len(sizes1) == 1:
            ax.bar(i - 0.2, sizes1[0] / 1e6, width=0.4, color='blue',edgecolor='black')  
        elif len(sizes1) == 2:
            ax.bar(i - 0.2, sizes1[0] / 1e6, width=0.4, color='blue',edgecolor='black') 
            ax.bar(i - 0.2, sizes1[1] / 1e6, width=0.4, color='blue', bottom=sizes1[0] / 1e6, hatch='//',edgecolor='black')  

        if len(sizes2) == 1:
            ax.bar(i + 0.2, sizes2[0] / 1e6, width=0.4, color='orange',edgecolor='black')  
        elif len(sizes2) == 2:
            ax.bar(i + 0.2, sizes2[0] / 1e6, width=0.4, color='orange',edgecolor='black')  
            ax.bar(i + 0.2, sizes2[1] / 1e6, width=0.4, color='orange', bottom=sizes2[0] / 1e6, hatch='//',edgecolor='black')

    # Set labels, title, and formatting
    ax.set_ylabel('Contig Size (Mb)')
    ax.set_xticks(range(len(contigs1)))
    ax.set_xticklabels(contigs1, rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(sys.argv[1], format="pdf")

# Create the combined barplot for both assemblies
create_combined_barplot(xie_et_al, isolate_2015T)
