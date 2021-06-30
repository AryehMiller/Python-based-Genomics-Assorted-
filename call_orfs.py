#!/usr/bin/env python3
"""Calling ORFs from a fasta file with multiple reading frames under a given condition

Usage: call_orfs.pyy <contigs fasta file>
Args: <contigs fasta file> a fasta file of sequences
Output: ORFs as amino acid residues and nucleotides in the format of FASTA files (all_orfs and all_proteins)

"""

## STANDARD IMPORT OF MODULES I USE

import sys
import csv
import re
import io
import os
import matplotlib.pyplot as plt
import numpy as np
import random
import string
import collections

## Quit if input parameters aren't correct

if len(sys.argv) != 2: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting a single fasta file"')

## READ IN FASTA
fasta_file = open(sys.argv[1], "r")

# Initializing contigs dictionary
fasta_dictionary = {}

# Initiate for loop to go through fasta_file
for line in fasta_file:
#Strip line
	line = line.rstrip()
#If statement-- denotes start of new sequence within the fasta file
	if line.startswith(">"):
#Assign the line to sequence name
		seq_name = line
#If line doesn't begin with >, then force all characters uppercase (precautionary) and put the sequence within the dictionary
	if not line.startswith(">"):
		line = line.upper()
	fasta_dictionary[seq_name] = line


## Create reverse complement dictionary for later use

reverse_dic = {'A':'T','T':'A','G':'C','C':'G'}
  
#Initiate dictionary
Combined_ORFs = {}
Combined_ORFs = collections.OrderedDict()

# Initiate for loop
## Includes code from http://doxey.uwaterloo.ca/tutorials/Python.html
## https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
for sequence_name,sequence_characters in fasta_dictionary.items():
	#Generate reverse complement
	contig_reverse_complement = "".join([reverse_dic[base] for base in reversed(sequence_characters)])
	#Generate list of sequence+reverse complement
	sequences = [sequence_characters, contig_reverse_complement]
	#Initiate a count for sequence orientation
	orientation = 0
	#Initiate for loop through all the sequences
	for sequence in sequences:
	#Add to the sequence count
		orientation += 1
		#Sliding window for reading frame
		for frame in range(0,3):
			#Assign true/false
			Containing_ORF = False
			#Initiate count as 0 to follow ORFs
			count = 0
			#Initiate for loop to move through sequence, 3 positions per time, dependent on reading frame
			for i in range(frame, len(sequence), 3):
			#Need to deposit the codon
				codon = sequence[i:i+3]
			#If start codon
				if codon == 'ATG':
			#And containing_ORF isn't true, then set to true and label the start position
					if Containing_ORF == False:
						Containing_ORF = True
						start_codon = i
			#If a stop codon appears, set the stop position to i and add to the count, the dictionary w/ the name+ the sequence
				if (codon in ("TAG","TAA","TGA")):
					if Containing_ORF == True:
						stop_codon = i
						count += 1
						Combined_ORFs[sequence_name+"_orientation_"+str(orientation)+"_rf_"+str(frame)+"_orf_"+str(count)] = sequence[start_codon:stop_codon+3]
						Containing_ORF = False

##Source of table: https://www.biostars.org/p/55851/
##Translate DNA sequence to amino acids
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
      'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
      'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W'}
      
      
## NOW STRINGENT CONDITIONS
ORFs = 0

## Prepare to write output files 
all_orfs = open("all_orfs.fna", "w")
all_proteins = open("all_proteins.faa", "w")

## Initiate for loop through all ORFs dictionary
for sequence_name,contig_orf in Combined_ORFs.items():
## If statement-- if the ORF is greater than or equal to 100 bp (not counting the stop codon)
	if len(contig_orf)-3 >= 100:
		# Add +1 to ORFs count
		ORFs += 1
		# Print to all_orfs
		print(sequence_name+'\n'+contig_orf, file = all_orfs)
		# Using genetic code at beginning of script, initialize protein sequence variable
		protein_sequence = ''
		#For loop to work through sequence, 3 positons at a time (not counting the stop)
		for i in range(0,len(contig_orf)-3,3):
		#If the ORF contig in the above sliding window matches the genetic code dictionary
			if contig_orf[i:i+3] in gencode:
				#Translate to amino acid and append
				protein_sequence += gencode[contig_orf[i:i+3]]
		#Print to all_proteins
		print(sequence_name+'\n'+protein_sequence, file = all_proteins)
	
## ORFs counted	
print("ORF count: " + str(ORFs))
## Close files
all_orfs.close()
all_proteins.close()
