#!/usr/bin/env python3
"""Calculates differences in regulatory potential for alternative-reference allele pairs using Log2FC scores

Usage: python3 count_barcodes.py <raw reads fastq> variant_to_barcode.txt > variant_fold_change.txt
Args: <raw reads fastq> fastq file
	 < variant to barcode file > a .txt file containing three columns with variant ID, reference barcodes (split by colons), and alternative barcodes (split by colons)
	 python3 analyze_MPRA.py pDNA_count.txt cDNA_count.txt filtered_variant_to_barcode.txt
Output: a table of variant, log2FC score, and p-value as determined by a Mann-U Whitney test
"""
import sys
import csv
import os
import math 
import scipy.stats

if len(sys.argv) != 4: #if the arguments isn't correct (4), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting three files-- two count files and variant-barcode file (in that order)')
	
## GENERATE NORMALIZED EXPRESSION VALUES FOR EACH BARCODE USING LOG2 RATIO OF CDNA/PDNA COUNTS
## norm expression = log2 (# of cDNA count_value / # pDNA count_value)

## INTIATE SEVERAL DICTIONARIES/LISTS/COUNTS FOR DOWNSTREAM USAGE

norm_expression_barcode = {}
cDNA_variant_reference_barcodes = {}
cDNA_variant_alternative_barcodes = {}
pDNA_variant_reference_barcodes = {}
pDNA_variant_alternative_barcodes = {}
pdna_dict = {}
cdna_dict = {}
variant_barcodes_reference = {}
variant_barcodes_alternative = {}
normalized_expression_reference = 0
norm_reference_count = 0
normalized_expression_alternative = 0
norm_alternative_count = 0
norm_values_reference = []
norm_values_alternative = []


## READ IN PDNA File
## Read in file-- standard line stripping and line splitting
with open(sys.argv[1]) as pDNA_file:
    for line in pDNA_file:
        line = line.rstrip("\n")
        #Define columns
        column = line.split("\t")
        #Get allele condition-- either ref or alt
        allele_condition = column[3]
        ## Extract variant identifier
        variant = column[2]
        ## Extract barcode and associated count
        barcode, count_value = column[0:2]
        pdna_dict[barcode] = count_value ## Add counts into dictionaries
        if allele_condition == "REF":
            if variant not in pDNA_variant_reference_barcodes:
                pDNA_variant_reference_barcodes[variant] = []
            pDNA_variant_reference_barcodes[variant].append(barcode)
        elif allele_condition == "ALT":
            if variant not in pDNA_variant_alternative_barcodes:
                pDNA_variant_alternative_barcodes[variant] = []
            pDNA_variant_alternative_barcodes[variant].append(barcode)

## READ IN CDNA File
## Read in file-- standard line stripping and line splitting
with open(sys.argv[2]) as cDNA_file:
    for line in cDNA_file:
        line = line.rstrip("\n")
        #Define columns
        column = line.split("\t")
        #Get allele condition-- either ref or alt
        allele_condition = column[3]
        ## Extract variant identifier
        variant = column[2]
        ## Extract barcode and associated count
        barcode, count_value = column[0:2]
        cdna_dict[barcode] = count_value ## Add counts into dictionaries
        if allele_condition == "REF":
            if variant not in cDNA_variant_reference_barcodes:
                cDNA_variant_reference_barcodes[variant] = []
            cDNA_variant_reference_barcodes[variant].append(barcode)
        elif allele_condition == "ALT":
            if variant not in cDNA_variant_alternative_barcodes:
                cDNA_variant_alternative_barcodes[variant] = []
            cDNA_variant_alternative_barcodes[variant].append(barcode)
            
## NORMALIZED EXPRESSION
## Compute the normalized expression value and add to dictionary "norm_expression_barcode"-- initiate for loop
for barcode in pdna_dict:
    if barcode in cdna_dict:
        normalized_expression = math.log2(float(cdna_dict[barcode])/float(pdna_dict[barcode]))
        norm_expression_barcode[barcode] = normalized_expression

## NOW NEED TO MAP BACK TO BARCODE FILE
## Open filtered variant-barcode file-- standard processing
with open(sys.argv[3]) as filtered_barcode_file:
    for line in filtered_barcode_file:
        line = line.rstrip("\n")
        column = line.split("\t")
        variants = column[0]
        reference_barcodes = column[1].split(":")
        alternative_barcodes = column[2].split(":")
        ## Extract ref and alt barcodes and add into respective dictionaries
        variant_barcodes_reference[variant] = reference_barcodes 
        variant_barcodes_alternative[variant] = alternative_barcodes

## NORMALIZED LOG FOLD CHANGE
## **** This is where I am having trouble- I'm not sure how to loop in each variant? *****
## Take norm_expression_barcode and assign back to ref or alt and identify fold change using log2FC
## Initiate for loop through variants
for variant in variant_barcodes_reference:
	for barcode in variant_barcodes_reference[variant]: ## Reference values for given variant
		normalized_expression_reference += norm_expression_barcode[barcode]
		norm_values_reference.append(norm_expression_barcode[barcode])
		norm_reference_count += 1 ## Add to count, if relevant
	for barcode in variant_barcodes_alternative[variant]: ## Alterative values for a given variant
		normalized_expression_alternative += norm_expression_barcode[barcode]
		norm_values_alternative.append(norm_expression_barcode[barcode])
		norm_alternative_count += 1 ## Add to count, if relevant
	mean_norm_expression_reference = normalized_expression_reference/norm_reference_count ## Average norm expression for reference
	mean_norm_expression_alternative = normalized_expression_alternative/norm_alternative_count ## Average norm expression for alternative
	log_trans_fold_change = mean_norm_expression_alternative - mean_norm_expression_reference ## Log2FC
	p_value = scipy.stats.mannwhitneyu(norm_values_alternative, norm_values_reference, alternative = "two-sided")[1] ## Determine significance using Mann-Whitney U test
	print(variant, log_trans_fold_change, p_value) ## Print out all three columns-- variant, the Log2FC value, and p-value from Mann-Whitney U test
