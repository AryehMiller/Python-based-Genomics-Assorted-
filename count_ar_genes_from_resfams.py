#!/usr/bin/env python3
"""Count genes in Resfams file

Usage: python3 count_ar_genes_from_blast.py <resfams .txt file>

Args: <resfams .txt file> a resfams output txt file
Output: the number of genes identified

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

## Initiate list
resfams = []

## Open the file
with open(sys.argv[1]) as resfams_file:
##For each line in the file
        for line in resfams_file:
#Split lines to read
            line = line.rstrip("\n")
#If the line doesn't start with "#"
            if line[0] != "#":
#Append the line the resfams
                resfams.append(line)
#Count the number of lines!           
print("Number of antibiotic resistant genes in Resfams: " + str(len(resfams)))
