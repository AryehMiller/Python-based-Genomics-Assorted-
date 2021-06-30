#!/usr/bin/env python3
"""Scan a DNA sequence to find putative binding sites

Usage: python3 demultiplex.py <scoring_matrix>

Args:
	<scoring_matrix> a csv file containing cell barcodes and cell tags
"""
import sys
import csv

# Check the correct number of command line arguments

if len(sys.argv) != 2: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting one file')

## Read in file 

csv_file = sys.argv[1]

#
##Control = ATGTTGC
## Treatment = GATTACA
## Non-determined = everything else
##If there is a one in the selected column, then add it to the count
##
#Source: https://stackoverflow.com/questions/9247241/python-algorithm-of-counting-occurrence-of-specific-word-in-csv
## How many Control cells

Control = csv.reader(open(csv_file)) #Open CSV
ctr = 0 #Start control count at 0
for record in Control: #for loop
    if record[1] == '1': #if observation = 1
        ctr += 1 #add 1 to the count
print("Control cells: " + str(ctr)) #print total

## How many Treatment cells
Treatment = csv.reader(open(csv_file)) #Open CSV
trt = 0 #Start treatment count at 0
for cell in Treatment: #for loop
    if cell[2] == '1': #if observation in column 3 = 1
        trt += 1 #add 1 to the count
print("Treatment cells: " + str(trt)) #print treatment count

## How many non-determined cells

ND = csv.reader(open(csv_file)) #Open CSV count
ND1 = 0 #Non-determined column 1 count at 0
ND2 = 0 #Non-determined column 2 count at 0
Zeros = 0 #These are the observations with zero
for cell in ND: #for loop
    if cell[3] == '1': #if observation in column 4 is "1"
        ND1 += 1 #add 1 to the non-determined count for ND1
    if cell[4] == '1': #if observation in column 5 is "1"
        ND2 += 1 #add 1 to the non-determined count for ND2
    if cell[1] == '0' and cell[2] == '0' and cell[3] == '0' and cell[4] == '0': #if all columns are 0, then add 1 to the zero count 
    	Zeros += 1 # Add 1 to zero count
ND_All = ND1+ND2+Zeros #Combine for non-determined total
print("Non-determined cells: " +str(ND_All)) #These are the total number of undetermined cells

## Edit distance of 1-- it would only affect the 3rd column with a 1-character distance for Control (column 4 has greater than 1-character diff)
Control = csv.reader(open(csv_file)) # Open CSV
ctr = 0 #Control starts at 0
EditDist = 0 #Edit distance group starts at 0
for record in Control: #For loop
    if record[1] == '1': #In column 2 if there's a 1
        ctr += 1 #+1 to the count
    if record[3] == '1': #if 1 in column 4
    	EditDist += 1 #Add 1 to the EditDist count
Edited_Dist = ctr + EditDist #Add together the control and the EditDist count to get the accommodated count since that's the only non-determined column eligible given an edit distance of 1
print("Control cells with edit distance of 1: " + str(Edited_Dist)) #Print it

print("Treatment cells with edit distance of 1: " + str(trt)) #Print it (will be same as original)

ND_All_Edit = ND_All - EditDist #The total non-determined with edit will be all - the edit distanced group
print("Non-determined cells with an edit distance of 1: " +str(ND_All_Edit)) #Hence we get the non-determined cells that fall within the edit threshold of 1



