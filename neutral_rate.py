#!/usr/bin/env python3

"""Calculates wobble proportion and conserved sequence ORFs

Usage: python3 neutral_rate.py <clustalw alignment file> <min bp conserved> <window length>
Args: <clustalw alignment file> = clustal alignment
 <min bp conserved> = an integer specifying the minimum number of conserved bp
 <window length> = a integer bp window length
Output: calculates fraction of wobble positions and promoter sequences within particular conservative criteria in an output file
"""
## Import relevant modules

import sys
import scipy.stats

## If number of input arguments isn't equal to 3, then exit

if len(sys.argv)!=4:
        sys.exit( __doc__)

## ARGUMENTS -- need min bp conserved threshold and the window length as integers
min_bp_conserved = int(sys.argv[2])
window_length = int(sys.argv[3])

## DOWNSTREAM USE IN EXTRACTING FASTA SEQUENCES
Skud_fasta = ""
Smik_fasta = ""
Scer_fasta = ""
Sbay_fasta = ""

## Open the alignment file and pull out each individual sequence and assign
with open(sys.argv[1]) as alignment:
	for line in alignment:
		if line.startswith("Skud"):
			Skud_fasta += line.strip().split()[1]
		if line.startswith("Smik"):
			Smik_fasta += line.strip().split()[1]
		if line.startswith("Scer"):
			Scer_fasta += line.strip().split()[1]
		if line.startswith("Sbay"):
			Sbay_fasta += line.strip().split()[1]

## CREATE DICTIONARY OF CODONS TO AMINO ACIDS
##Source of table: https://www.biostars.org/p/55851/
## NOTE -- need to use "*" for stops since gaps (-) are present in alignment 

gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}


## MODIFIED FROM CALL_ORFS.PY IN ASSIGNMENT 11
## Includes code from http://doxey.uwaterloo.ca/tutorials/Python.html
## https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python					
						
## LOOP THROUGH ORFS IN 3 BP WINDOW TO FIND CONSERVED REGIONS
## SAME GENERAL STRUCTURE AS CALL_ORFS.PY IN ASSIGNMENT 10
Total_sequence_len = 0 ## Need further down for count
for reading_frame in range(0,3): ## 3 bp window
	relevant_ORF = False ## Initialize false for ORF
	for j in range(reading_frame,len(Smik_fasta),3): ## Loop through reading frame for length of  individual in fasta
		## Initiate if statement for each AA window in the fasta file sequences
		if (Skud_fasta[j:j+3] in gencode) and (Smik_fasta[j:j+3] in gencode) and (Scer_fasta[j:j+3] in gencode) and (Sbay_fasta[j:j+3] in gencode):
			Skud_AA = gencode[Skud_fasta[j:j+3]] ## Translate to AA
			Smik_AA = gencode[Smik_fasta[j:j+3]] ## Translate to AA
			Scer_AA = gencode[Scer_fasta[j:j+3]] ## Translate to AA
			Sbay_AA = gencode[Sbay_fasta[j:j+3]] ## Translate to AA
			if 'M' == Skud_AA == Smik_AA == Scer_AA == Sbay_AA: ## if shared start
				if relevant_ORF == False:
					relevant_ORF = True
					pos_start = j ## Need to keep track
			if '*' == Skud_AA == Smik_AA == Scer_AA == Sbay_AA: ## If there's a stop (*) 
				if relevant_ORF == True:
					pos_stop = (j + 3) ## Assign to stop pos
					if (pos_stop - pos_start) > Total_sequence_len: ## If the stop position - start position is greater than the total length of the sequence
						Total_sequence_len = (pos_stop - pos_start)
						START_ORF = pos_start
						STOP_ORF = pos_stop
					relevant_ORF = False

## IDENTIFY PROPORTION OF CONSERVED WOBBLE POSITIONS
## SAME GENERAL STRUCTURE AS CALL_ORFS.PY IN ASSIGNMENT 10
wobbles_count = 0 ## Initiate count for wobbles_count at 0
conserved_count = 0 ## Initiate count for conserved positions at 0
for j in range(START_ORF, STOP_ORF, 3): 
## Initiate for loop to go through each fasta and add to either wobble or conserved positions count
## Need to deposit the codons
	Skud_AA = Skud_fasta[j:j+3]
	Smik_AA = Smik_fasta[j:j+3]
	Scer_AA = Scer_fasta[j:j+3]
	Sbay_AA = Sbay_fasta[j:j+3]
	if gencode[Skud_AA] == gencode[Smik_AA] == gencode[Scer_AA] == gencode[Sbay_AA] and gencode[Skud_AA] != 'M' and gencode[Skud_AA] != 'W':
		if Skud_AA[0:2] == Smik_AA[0:2] == Scer_AA[0:2] == Sbay_AA[0:2]:
			wobbles_count += 1
			if Skud_AA == Smik_AA == Scer_AA == Sbay_AA:
				conserved_count += 1

print('Proportion conserved: ' + str(conserved_count/wobbles_count))

## BINOMIAL TEST
## Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom_test.html
print("Binomial test: ") 
print(scipy.stats.binom_test(8, n=10, p=0.44, alternative='greater')) ## 8 provides a p-value less than 0.05

## IDENTIFY PROMOTER SEQS
## SAME GENERAL STRUCTURE/LOGIC AS CALL_ORFS.PY IN ASSIGNMENT 10
list_of_promoters = [] ## Initiate list of prom positions
for prom in range(0, START_ORF - window_length): ## Initiate for loop 
	number_of_conserved_bases = 0 ## Initiate count of the conserved bases at 0
	for bp in range(prom, prom + window_length): ## Initiate for loop w/ respect to promoters
		if Skud_fasta[bp] == Smik_fasta[bp] == Scer_fasta[bp] == Sbay_fasta[bp]: ## If shared across sequence bp, add one to the conserved base count
			number_of_conserved_bases += 1 ## Add 1
	if number_of_conserved_bases >= min_bp_conserved: ## Conditional-- if the number of conserved bases is greater than or equal to the specified minimum base threshold, then append to the list
		list_of_promoters.append(prom)

## WRITE OUTPUT FILE WITH START, STOP, AND SEQUENCE NAMES
## Unsure of how to exactly go about doing this-- need to keep track of start and stop positions, with regard to sequence...
with open('S_cer_conserved.txt', 'w') as text_file:
	print("START", "STOP", "SEQUENCE", file = text_file)
