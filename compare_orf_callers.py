#!/usr/bin/env python3
"""Compare ORF callers

Usage: python3 compare_orf_callers.py <fasta 1> <fasta 2>
Args: <fasta 1> a standard fasta file of ORFs
	  <fasta 2> a standard fasta file of ORFs
Output: Several counts about the ORFs shared and unique to the two files

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

if len(sys.argv) != 3: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting two fasta files')

## FIRST STEP: OPEN UP call_orfs.py OUTPUT files and the metamark files (here I'm using the protein sequences)

#Initialize lists to hold
all_proteins = []
mgm_orfs = []

## ALL_PROTEINS-- open up with standard previously used method-- append characters to all_proteins
with open(sys.argv[1]) as fasta_file:
        for line in fasta_file:
            line = line.rstrip("\n")
            if line[0] != ">":
                all_proteins.append(line)

## MetaMarkGene-- open up with a bit of a different approach but similar nonetheless-- see https://www.biostars.org/p/18129/
with open(sys.argv[2]) as fasta_file:
    ORFs = ""
    for line in fasta_file:
        line = line.rstrip("\n")
        if line[0] != ">":
            ORFs += line
        else:
            if ORFs != "":
                mgm_orfs.append(ORFs)
                ORFs = ""
    mgm_orfs.append(ORFs)

## COMPARE AND FIND MATCHES
#https://stackoverflow.com/questions/13323851/python-3-counting-matches-in-two-lists-including-duplicates/13323960
#https://stackoverflow.com/questions/35713093/how-can-i-compare-two-lists-in-python-and-return-not-matches
#SHARED
print("Number of ORFs in mgm_orfs.faa: " + str(len(mgm_orfs))) ## Print # of ORFs in mgm_orfs
print("Number of ORFs in all_proteins.faa: " + str(len(all_proteins))) ## Print # of ORFs in all_proteins
print("Number of shared ORFs: " + str(len([w for w in mgm_orfs if w in all_proteins]))) ## This pulls the shared ORFs among the fastas
print("Number of ORFs unique to mgm_orfs.faa: " + str(len([w for w in mgm_orfs if w not in all_proteins]))) ## This pulls the # of ORFs unique to the mgm
print("Number of ORFs unique to all_proteins.faa: " + str(len([w for w in all_proteins if w not in mgm_orfs]))) ## This pulls the # of ORFs unique to all_proteins
