#!/usr/bin/env python3

"""
usage: python3 analyze_WGBS_methylation.py <WGBS bed>
Arguments: <WGBS bed> a WGBS bed file 
Output: 4-column bedfile: 1) chromosome, 2) start, 3) stop, and 4) methylation level.
"""
#Import modules

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
BGM_WGBS = open(sys.argv[1],"r") 
filepath = sys.argv[1]

#Get the name
name = os.path.basename(filepath)
#Get the pruned name (w/o file ending)
short_name = os.path.splitext(name)[0]
##Define the dynamic output file name in each case
out_name = short_name + "_CpG_methylation.bed"
#Object to write
out_object = open(out_name, "w")

#Store values > 0
CPG_METHYL_LIST = []
#Store the coverage list
store_CpG_coverage = []
#Begin the count at zero
CPG_METHYL_LIST_count_zero = 0
#Begin the count at zero
island_count = 0
#Initiate for loop to run through the .bed file
for line in BGM_WGBS: #for each line in the .bed file, 
	each_line = line.strip().split() #strip white space
	CPG_COVERAGE = (float(each_line[3])+float(each_line[4])) #grab values from the two relevant columns in the .bed file and put into CpG_coverage using a float operation
	store_CpG_coverage.append(CPG_COVERAGE) #Append to cpg coverage list

	if each_line[3] == '0' and each_line[4] == '0': #if and statement to select all values (essentially filter out zeros) values that aren't equal to zero in both relevant columns-- this is kind of a workaround to doing something otherwise in pandas
		CPG_METHYL_LIST_count_zero += 1 #add to the count
	else: #otherwise
		island_count += 1 #add to the island count
		island_methylation = (float(each_line[3])/(float(each_line[3])+float(each_line[4]))) #Here's where each actually implement the referenced equation in the lab instructions-- C bases/C+T bases-- to calculate the CpG methylation level
		CPG_METHYL_LIST.append(island_methylation) #append
		print(str(each_line[0]) + '\t' + str(each_line[1]) + '\t' + str(each_line[2]) + '\t' +  str(island_methylation), file = out_object) # print out new table using str to make it look nice
out_object.close() #close the file


## Plots the distribution with a histogram of CpG methylation levels as <WGBS bed basename>_methylation_distribution.png.
## Using this source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html
plt.hist(CPG_METHYL_LIST, bins = 100, color = "lightskyblue",edgecolor='black', linewidth=0.5) ## using a histogram and binwidth at 100-- note, however, that one could also make this into a density distribution for a different visualizaiton
plt.title('CpG Methylation Distribution in ' + short_name) #Define main title
plt.xlabel("CpG Methylation Value") #Define x axis title
plt.ylabel("Site Frequency") #Define y axis title
plt.show() #show the plot
plt.savefig(short_name+'_methylation_distribution.png') #Save the plot

#Plots the distribution of read coverage for all CpGs for coverages between 0X and 100X as <WGBS bed basename>_CpG_coverage_distribution.png.
plt.hist(store_CpG_coverage, bins = 5000, color = "lightskyblue",edgecolor='black', linewidth=0.5)  ## using a histogram and bin # at 5000-- note, however, that one could also make this into a density distribution for a different visualizaiton
plt.title('Distribution of CpG Coverage beteen 0x and 100x in ' + short_name) #Define main title
plt.xlabel("CpG Methylation Coverage") #Define x axis title
plt.ylabel("Frequency") #y axis title 
plt.xlim([0, 1000]) #x lim
plt.show() #show the plot
plt.savefig(short_name+'_CpG_coverage_distribution.png') #Save the plot

##Calculate 0x coverage fractions
print('the 0x coverage fraction is:' +str(CPG_METHYL_LIST_count_zero/island_count))#Almost 1%

