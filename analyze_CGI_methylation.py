#!/usr/bin/env python3

"""
usage: python3 analyze_CGI_methylation.py <average CGI methylation bed>
Arguments: <WGBS bed> a BED file with the columns 1) chromosome, 2) start, 3) stop, 4) CGI name, and 5) average CpG methylation in CGI
Output: plot of methylation distribution
"""
import sys
import os
import matplotlib.pyplot as plt

##Check input parameters as two

if len(sys.argv) != 2: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting two')

### Make dynamic naming system for files

### Saves the input arguments as variables

### Open the first argument file and get the path
Mean_cpg = open(sys.argv[1],"r")
filepath = sys.argv[1]

#Get the name
name = os.path.basename(filepath)
#Get the pruned name (w/o file ending)
shortname = os.path.splitext(name)[0]
#Initiate empty list to hold values
Mean_CpG_Value = []

#For loop to obtain 5th column using a for loop
for line in Mean_cpg:
	each_line = line.strip().split() #strip white space
	Mean_CpG_Value.append(float(each_line[5])) #extract and append 5th column values

#Histogram
## Using this source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html
plt.hist(Mean_CpG_Value, bins = 100, color = "lightskyblue",edgecolor='black', linewidth=0.5) #Using a histogram w/ 100 bins
plt.title('Methlylation distribution for \n' + shortname + '_CpG_methylation.bed') #Define main title-- use \n which places the following line below (useful for saving space)
plt.xlabel("Mean CpG Methylation Value") #x axis label
plt.ylabel("Frequency") #y axis label
plt.show() #show plot
plt.savefig(shortname+'_distribution.png') #Save the plot

