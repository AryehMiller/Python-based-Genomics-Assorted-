#!/usr/bin/env python3
"""Interrogate genomic structural variation from a VCF file

Usage: python3 count_gv.py <SNV_indel VCF> <SV VCF> 

Args:
	<SNV_indel VCF> a snv VCF file
	<SV VCF> a structural variants VCF file
	
"""
import sys
import csv
import re
import io
import os
import matplotlib.pyplot as plt


# Check the correct number of command line arguments
if len(sys.argv) != 3: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting two files')

## Read in file 

SNV_indel = sys.argv[1]
SV_VCF = sys.argv[2]

#Set up counts 
deletions = 0
insertions = 0
duplications = 0
inversions = 0
copy_number_variations = 0
tandem_duplications = 0
mobile_element_deletion = 0
mobile_element_insertion = 0
breakends = 0

## COUNT AND PRINT SNV INDELS AND TOTAL SNVS

#Set up the variant lists
indel_length = []
mobile_element_insertion_length = []
deletion_length = []
cnv_length = []

#Set up Indel/SNV counts
indel_count = 0
snv_count = 0

#count SNVs and indels
with open(SNV_indel) as SNV_indel_file: #Open the first file
        for line in SNV_indel_file: #Go through each line
            line = line.rstrip("\n")
            columns = line.split("\t")
            if line[0] != "#": #ignore first
                na12878 = columns[10] #Extract na12878
                if "0/0" not in na12878 and "./." not in na12878: #If 0/0 or missing ./. isn't present
                    if len(columns[3]) == 1 and len(columns[4]) == 1: #SNVs
                        snv_count += 1
                    else: #Indels
                        indel_count += 1 #Add 1 to indel count
                        indel_length.append(abs(len(columns[3]) - len(columns[4]))) #Frequency for plotting
                        
print("Indel count: " + str(indel_count)) ## PRINT total indels 
print("Total SNV count: " + str(snv_count)) ## PRINT total SNVs

## COUNT STRUCTURAL VARIANTS

with open(SV_VCF) as StructVarFile: #Open up 2nd file
        for line in StructVarFile: #Go through each line
            line = line.rstrip("\n")
            columns = line.split("\t")
            if line[0] != "#": #Ignore first
                na12878 = columns[10] #Extract na12878 column
                if "0/0" not in na12878 and "./." not in na12878: #If 0/0 or missing ./. isn't present
                    information =  columns[7] #Extract 8th column
                    information_split = information.split(";") #Split at semi-colon
                    for i in information_split: #For loop in information_split
                    
                    ## EXTRACT STRUCTURAL VARIANTS PER https://samtools.github.io/hts-specs/VCFv4.2.pdf
                    ## Add to count if structural variant type is present
                        if "SVLEN" in i:
                            svlen = abs(int(i.replace("SVLEN=", ""))) #absolute value of length
                            break
                    #BND
                    if "SVTYPE=BND" in information:
                        breakends += 1
                        
                    #Tandem duplications
                    elif "SVTYPE=DUP:TANDEM" in information:
                        tandem_duplications += 1   
                        
					#Mobile element insertion
                    elif "SVTYPE=INS:ME" in information:
                        mobile_element_deletion += 1 
                        
                    #Mobile element deletion
                    elif "SVTYPE=INS:ME" in information:
                        mobile_element_deletion += 1 
                        
                    #MEI
                    elif "SVTYPE=MEI" in information:
                        mobile_element_insertion += 1
                        mobile_element_insertion_length.append(svlen)
                    #DEL
                    elif "SVTYPE=DEL" in information:
                        deletions += 1
                        deletion_length.append(svlen)
                    #DUP
                    elif "SVTYPE=DUP" in information:
                        duplications += 1
                    #INV
                    elif "SVTYPE=INV" in information:
                        inversions += 1
                    #CNVs
                    elif "SVTYPE=CNV" in information:
                        copy_number_variations += 1   

#Print the number for each category
print("Deletions count: " + str(deletions))
print("Insertions count:" + str(insertions))
print("Duplications count: " + str(duplications))
print("Inversions count: " + str(inversions))
print("Copy number variants: " + str(copy_number_variations))
print("Tandem duplications: " + str(tandem_duplications))
print("Mobile element deletions: " + str(mobile_element_deletion))
print("Mobile element insertions: " + str(mobile_element_insertion))

## Proportion of genomic variants that are SVs proportion of genomic variants that are SVs
Total_SVs = deletions +insertions +duplications +inversions +copy_number_variations +tandem_duplications +mobile_element_deletion +mobile_element_insertion +breakends ## Combine variants for total SVs
Total_Genomic_Variants = snv_count + indel_count + deletions +insertions +duplications +inversions +copy_number_variations +tandem_duplications +mobile_element_deletion +mobile_element_insertion +breakends ## Total genomic variants including SNVs
Proportion_SVs = Total_SVs/Total_Genomic_Variants*100 #What proportion are major structural variants-- I want a percentage, so I multiply by 100
print("Proportion of genomic variants that are structural variants: " + str(round(Proportion_SVs,3)) + "%") #Print and round

## HISTOGRAMS

## Indels
## Using this source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html
plt.hist(indel_length, bins = 500, color = "lightskyblue",edgecolor='black', linewidth=0.5) ## using a histogram and binwidth at 500-- note, however, that one could also make this into a density distribution for a different visualizaiton
plt.title('Distribution of small indels in NA12878') #Define main title
plt.xlabel("Length") #Define x axis title
plt.ylabel("Number of Indels") #Define y axis title
plt.xscale('log') #Log scale
plt.show() #show the plot
plt.savefig('histogram_indels.png') #Save the plot

## MEIs
## Using this source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html
plt.hist(mobile_element_insertion_length, bins = 500, color = "lightskyblue",edgecolor='black', linewidth=0.5) ## using a histogram and binwidth at 500-- note, however, that one could also make this into a density distribution for a different visualizaiton
plt.title('Distribution of MEIs in NA12878') #Define main title
plt.xlabel("Length") #Define x axis title
plt.ylabel("Number of MEIs") #Define y axis title
plt.xscale('log') #Log scale
plt.show() #show the plot
plt.savefig('histogram_meis.png') #Save the plot


## Deletions
## Using this source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html
plt.hist(deletion_length, bins = 500, color = "lightskyblue",edgecolor='black', linewidth=0.5) ## using a histogram and binwidth at 500-- note, however, that one could also make this into a density distribution for a different visualizaiton
plt.title('Distribution of deletions in NA12878') #Define main title
plt.xlabel("Length") #Define x axis title
plt.ylabel("Number of Deletions") #Define y axis title
plt.xscale('log') #Log scale
plt.show() #show the plot
plt.savefig('histogram_deletions.png') #Save the plot


## IF one wanted to use grep for some other analysis
#grep -i "NA12878" /home/assignments/assignment8/sv.reclassed.filtered.vcf > NA12878.sv.reclassed.filtered.vcf