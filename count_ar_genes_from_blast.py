#!/usr/bin/env python3
"""Examine BLAST output

Usage: python3 count_ar_genes_from_blast.py <blast output>
Args: <blast output> a textfile in tabular format from BLAST
Output: the number of hits >80% AA identity and > 85% coverage of subject sequence
"""
## STANDARD IMPORT OF MODULES I USE

import sys
import csv
import re
import io
import os
import matplotlib.pyplot as plt
import numpy as np
import random
import string
import collections


if len(sys.argv) != 2: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting one input file')


## APPLY FILTERS

## Open and read file
blast_file = open(sys.argv[1], "r")

## Initiate gene count
gene_count = 0
## Initiate for loop through BLAST file
for line in blast_file:
#Strip line
	line = line.rstrip()
#Pull specific columns
	seq_coverage =  float(line.split()[12]) ##Sequence coverage (need float)
	AA_identity = float(line.split()[2]) ## Amino acid identity (need float)
	sequence_length = float(line.split()[3]) ##Length of sequence (need float)
	## If the amino acid identity is > 80% and the coverage is greater than a proportion of 0.85 (seq length / seq coverage)
	if AA_identity > 80 and sequence_length/seq_coverage > 0.85:
	#Add one to the count
		gene_count+=1

print("Sequences passing both filters: " + str(gene_count))

