#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 19:03:17 2025

@author: yomna
"""


import pysam
import sys
bam_file = sys.argv[1]
chromosome= sys.argv[2]
start_pos = int(sys.argv[3])
end_pos = int(sys.argv[4])
samfile = pysam.AlignmentFile(bam_file, "rb")
longest_read = []
longest_length = 0
for read in samfile.fetch(chromosome, start_pos, end_pos):
    if read.query_length > longest_length: #and (read_start >= start_pos and read_end <= end_pos):
        longest_length = read.query_length
        longest_read = read
samfile.close()

print("longest_read", longest_read.query_name)