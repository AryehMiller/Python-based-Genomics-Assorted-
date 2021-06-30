#!/usr/bin/env python3
"""
Usage: python3 generate_promoters.py <standard bed file>
Arguments: <standard bed file>  = standard bed file containing labels and gene coordinates
Output: 5-column bedfile: 1) chromosome, 2) start, 3) stop, 4) gene name, and 5) strand
"""
# Import modules needed for program
import sys
import os
import matplotlib.pyplot as plt


if len(sys.argv) != 2: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting two')


#Read in bed file
Promoters = open(sys.argv[1],"r")
filepath = sys.argv[1]

# Creating dynamic output filenames
#Get the name
name = os.path.basename(filepath)
#Get the pruned name (w/o file ending)
shortname = os.path.splitext(name)[0]
##Define the dynamic output file name in each case
out_name = shortname + "_promoters.bed"
#Object to write
out_object = open(out_name, "w")

#Initiate For loop for each line in bed file
for line in Promoters:
	row = line.strip().split() # Strip whitespace
	START = int(row[1]) - 850 #start position by subtracting 850 from initial position
	#stop position by subtracting 150 from initial position
	STOP = int(row[1]) - 150 #Stop position
	print(str(row[0]) + '\t' + str(START) + '\t' + str(STOP) + '\t' +  str(row[3]) + '\t' + str(row[5]), file = out_object)# Print to the output file
out_object.close()#Close object