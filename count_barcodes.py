#!/usr/bin/env python3
"""Filters out reads that don’t perfectly match a barcode, and counts the total number of perfectly matching reads for each barcode

Usage: python3 count_barcodes.py <raw reads fastq> variant_to_barcode.txt
Args: <raw reads fastq> fastq file
Output: # of filtered read, and a file with barcode-variant reads counts
"""
import sys
import os
import csv

if len(sys.argv) != 3: #if the arguments isn't correct (3), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting two input files-- a .fq file and a variant to barcode .txt file')
	
## READ IN RELEVANT FILES

fastq_file = sys.argv[1]
variant_to_barcode = sys.argv[2]

### MAKE DICTIONARIES FOR REFERENCE AND ALT BARCODES WITH ASSOCIATED IDS

## Initiate several empty dictionaries for barcodes
reference_barcode_dict = {}
alt_barcode_dict = {}
barcode_total_count_dict = {}
filt_read_count = 0 ## Will use this downstream for filtered reads in writing new file

## PARSE VARINT-BARCODE FILE
## Same general structure as previous script (filter_variants.py), except at the end we're storing the respective reference and alternative barcodes in the above dictionaries
with open(variant_to_barcode) as variant_to_barcode_file:
    for line in  variant_to_barcode_file:
        line = line.rstrip("\n")
        column = line.split("\t")
        variant = column[0]
        reference_barcode = column[1].split(":")
        alternative_barcode = column[2].split(":")
        reference_barcode_dict[variant] = reference_barcode
        alt_barcode_dict[variant] = alternative_barcode

## PARSE FASTQ FILE

with open(fastq_file) as fq_file:
	#note-- we'll need to keep track of the line number, hence line_number count
    line_number = 0
    for line in fq_file: 
        line_number += 1
        line = line.rstrip("\n")
        if line_number % 2 == 0 and line_number % 4 != 0: ## Within the .fq file, examine the sequence line only
            barcode = line[14:23] ## Extract 9 bp barcode
            if barcode not in barcode_total_count_dict: ## Check if it's already within the dictionary, add if not
                barcode_total_count_dict[barcode] = 0
            barcode_total_count_dict[barcode] += 1

## TIME TO WRITE A NOVEL FILE-- this is why I was forced to copy the .fq files to my own working directory because of permissions issues in the master directory for the whole class

## Prepare file naming system using .replace-- we want to have fastq file name + _count.txt
count_file = fastq_file.replace(".fq", "") + "_count.txt"

## MATCH UP DICTIONARY VALUES/KEYS WITH FASTQ FILES AND KEEP RUNNING COUNT

## Prepare to write a novel file for each of the fastq files, respectively.
with open(count_file, "w") as novel_count_file: ## Open count file
    for barcode in barcode_total_count_dict: #initiate for loop-- for barcode in the total dictionary
        allele_condition = "NA"
        condition = "NA"
        for variant in reference_barcode_dict:
            if barcode in reference_barcode_dict[variant]: ## if barcode matches filtered version (REF)
                allele_condition = "REF"
                condition = variant
        if allele_condition == "NA": ## if not reference
            for variant in alt_barcode_dict:
                if barcode in alt_barcode_dict[variant]: ## if barcode matches filtered version (ALT)
                    allele_condition = "ALT"
                    condition = variant
        if allele_condition != "NA": ## If variant condition not found
        	## print non-filtered barcodes
            print(barcode, barcode_total_count_dict[barcode], condition, allele_condition, sep = "\t", file = novel_count_file)
        else:
        	## Add in filtered reads
            filt_read_count += barcode_total_count_dict[barcode]
##Print number of filtered out reads-- we can just calculate this using  a simple awk command otherwise (see my README doc)
print("Filtered out reads: ", filt_read_count)

