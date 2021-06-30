#!/usr/bin/env python3

"""Filters and calculates variant score means

Usage: python3 topscores.py pten_mave_scores.csv pten_variants.txt
Args: <pten_mave_scores.csv> a csv file containing variant names and corresponding scores in columns 3 and 4
	  <pten_variants.txt> a text file containing a single column of variant names of interest
Output: Mean scores and corresponding scores for list of candidate variants
"""
## Import relevant modules

import sys
import csv
from itertools import islice
from statistics import mean

## If number of input files isn't equal to 2, then exit

if len(sys.argv)!=3:
        sys.exit( __doc__)
    

def load_mave_data(mave_file):
    ### Returns a dictionary of mave data variant names and scores
    ### Read in variant names and scores from mave file
    ### Remember to skip over the first 5 rows
    ### Store the data from the 3rd and 4th columns in the dict 'scores'
    ### REMEMBER to convert scores from strings to floats!
    # define dictionary
    # skip through first rows:
	
    # read in file data:    
	numbers = [] ## Initiate list for 4th column
	name = [] ## Initiate list for 3rd column
	scores = {} ## define 'scores' dictionary
	
	## Open MAVE file and isolate scores column
	with open(sys.argv[1]) as f:
	    reader = csv.reader(f) ## Open the CSV file and assign to 'reader'
	    reader = islice(reader, 4, None) ## Read in everything after 4th row
	    for row in reader: ## Initiate for loop-- extract 4th column and store in numbers
	    	numbers = [i[3] for i in reader]
	## Open MAVE file and isolate variant name column-- same code, but for the 3rd column (variant names) and store in 'name'
	with open(sys.argv[1]) as f:
	    reader = csv.reader(f) 
	    reader = islice(reader, 4, None)
	    for row in reader:
	    	name = [i[2] for i in reader]
	
	## Convert score strings into float values    	
	score_values = [float(i) for i in numbers]
	scores = dict(zip(name, score_values)) ## Assemble and return dictionary 'scores' with names and score values 
	return scores
    
def score_variants(variant_file, scores):
    ### Function should return a *sorted* list of variants and scores
    ### AND the mean of the variant scores
    ### Read in the variant file as a list
    ### Look up each variant score in the dictionary 'scores'
    ### Save variant/score pairs as a *list of tuples* (variant, score)
    ### Sort the list of tuples by score
    
    # define variant list
	dictionary_all_variants = [] 
	dictionary_all_variants = [(k, v) for k, v in scores.items()] ## Assemble dictionary for all variants
	sorted_variants_ALL = sorted(dictionary_all_variants, key = lambda x: x[1]) ## Sort dictionary of all variants
	variants = [] ## Make empty list 'variants'
	mave_mean = [] ## Make empty list 'mave_mean'
	variants_mean = [] ## Make empty lsit 'variants_mean'
	## Get list of 100 important candidate variants
	with open(sys.argv[2]) as variant_list: ## Open up the variant list .txt file
		for line in variant_list: ## For each line in the file
			variants.append([ x for x in line.split()]) ## Append to 'variants'
			column1 = [ x[0] for x in variants] ## Assign such to "column1"
    
    ##MAKE FILTER
	filter = column1 ## Assign the 100 variants to 'filter'
	filter_set = set(filter) ## set filter
	tuples_filtered_one_hundred = [tup for tup in sorted_variants_ALL if tup[0] in filter_set] ##Create tuple of one hundred variants based on filter
	
    # sort the list of tuples by score (using a lambda function)
	sorted_variants = sorted(tuples_filtered_one_hundred, key = lambda x: x[1])
    
    # calculate the mean of the variant scores:
    ## Mean of variant data
	variants_mean = mean(value[1] for value in sorted_variants)
	variant_file.close()
	return(variants_mean, sorted_variants)

## Main script
mave_file = open(sys.argv[1])
variant_file = open(sys.argv[2])

scores = load_mave_data(mave_file)
variants_mean, sorted_variants = score_variants(variant_file, scores)

### Calculate the mean of the MAVE data:
mave_mean = [] ## Empty list
dictionary_all_variants = [(k, v) for k, v in scores.items()] ## Dictionary of all variants with scores
mave_mean_dataset = sorted(dictionary_all_variants, key = lambda x: x[1]) ## Sort that dataset (not necessary, though)
mave_mean = mean(value[1] for value in mave_mean_dataset) ## Calculate mean of MAVE scores
print('Mave mean', mave_mean, 'Variants mean', variants_mean, '\n') ## Print MAVE mean and variants mean
## Print 100 variant scores
for name, score in sorted_variants:
    print(name, score, '\n')    





