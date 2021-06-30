#!/usr/bin/env python3
"""Filter variants assigned that pass a 7-barcode filter

Usage: python3 filter_variants.py <variant-barcode file>  > filtered_variant_to_barcode.txt
Args: <variant-barcode file> a file comprising variant IDS in one column, followed by REF and ALT barcodes in the next two. Each barcode is separated by a colon.
Output: a filtered variant-barcode file according to threshold (here it is 7 or greater barcodes)
"""
import sys
import csv

# Check the correct number of command line arguments

if len(sys.argv) != 2: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting one file')

## READ IN .TXT FILE
variants_file = open(sys.argv[1])

### Filter based on colons-- if line has at greater than or to 7 ":" (denoting barcodes), then include
with variants_file as variants_file_file: #Open the first file
		#next(variants_file_file) ## this will exclude the headers in entirety
		header = True
		## Just skip the header
		for line in variants_file_file: ## Go through each line
			line = line.rstrip("\n") ## Line stripping
			## we need to take care of the header issue-- want to exclude them in our iteration through barcodes
			if header:
				print(line)
				header = False
				next
			## Split table up in columns
			columns = line.split("\t")
			## Shaft columns 1 and 2 into list, split by colons
			REF_barcode = columns[1].split(":") ## Reference barcode column
			ALT_barcode = columns[2].split(":") ## Alternative barcode column
			## If there are 7 or more barcodes, then print that given variant
			if(len(REF_barcode) >= 7 and len(ALT_barcode) >= 7):
				print(line)
        
   