#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 19:00:11 2024

@author: yomna
"""
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.cm import get_cmap
import sys
# Step 1: Parse the PAF file to get read segments
def plot_rDNA(paf_file,output=None):
    read_segments = {}
    reference_names = []  # We'll extract these dynamically from the PAF file
    with open(paf_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            read_name = fields[0]
            read_start = int(fields[2])
            read_end = int(fields[3])
            strand= fields[4]
            ref_name = fields[5]
            if ref_name not in reference_names:
                reference_names.append(ref_name)
            
            # Store the segment information
            if read_name not in read_segments:
                read_segments[read_name] = []
            read_segments[read_name].append({
                "start": read_start,
                "end": read_end,
                "strand":strand,
                "ref_name": ref_name
            })
    for read_name in read_segments:
        read_segments[read_name].sort(key=lambda x: x["start"])  # Sort by "start" value    
    reference_file = sys.argv[1]
    reference_sequences_1 = list(SeqIO.parse(reference_file, "fasta"))
    num_regions = len(reference_names)
    cmap = get_cmap("tab20", num_regions)  # Use a qualitative colormap
    region_colors = {ref_name: cmap(i) for i, ref_name in enumerate(reference_names)}
    
    region_descriptions = {ref.id: ref.description for ref in reference_sequences_1}
    region_colors_1 = {ref_name: region_descriptions.get(ref_name, "No description available") for ref_name in reference_names}
    plt.figure(figsize=(30, len(read_segments) * 0.5 + 3))  # Increase figure height dynamically
    y_pos = 0
    for read_name, segments in read_segments.items():
        y_pos += 1
        for segment in segments:
            # Plot the rectangle (bar)
            plt.barh(
                y_pos,
                width=segment["end"] - segment["start"],
                left=segment["start"],
                color=region_colors[segment["ref_name"]]
                #edgecolor="black"
            )
    
            # Add strand label (+ or -) above each block
            plt.text(
                (segment["start"] + segment["end"]) / 2,  # Position text at the center of the rectangle
                y_pos + 0.5,  # Move label slightly above the rectangle
                segment["strand"],  # Label showing the strand
                ha="center",  # Center align the text
                va="center",  # Center vertically the label
                fontsize=10,  # Font size for the label
                color="black"
            )
    
    plt.xlabel("Position on Read")
    plt.ylabel("Reads")
    plt.yticks([])  
    plt.tight_layout(rect=[0, 0.1, 1, 1])     
    handles = [plt.Line2D([0], [0], color=color, lw=4) for color in region_colors.values()]
    plt.legend(
        handles,
        region_colors_1.values(),
        loc="upper center",
        bbox_to_anchor=(0.5, -0.2),  # Moves legend completely below the plot
        ncol=5,
        title="Reference Regions"
    )
    
    plt.savefig(output,format="pdf")
    
paf_file = sys.argv[2]
plot_rDNA(paf_file,output=sys.argv[3])
