#!/usr/bin/env python3
"""Quantify (count) genotypes in a VCF file

Usage: python3 quantify_genotype.py <SNV_indel VCF>

Args:
	<SNV_indel VCF> a VCF file for which one is interested in generating SNV/Indel genotype counts
"""
import sys
import csv
import re
import io
import os
import matplotlib.pyplot as plt


# Check the correct number of command line arguments

if len(sys.argv) != 2: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting one file')

## Read in file 

SNV_indel = sys.argv[1]

#Set up counts 

#Set up counts 
SNV_homozygous_ref = 0
SNV_heterozygous = 0
SNV_homozygous_alternative = 0
SNV_missing_genotypes = 0
INDELS_homozygous_ref = 0
INDELS_heterozygous = 0
INDELS_homozygous_alternative = 0
INDELS_missing_genotypes = 0
indel_count = 0
snv_count = 0

with open(SNV_indel) as SNV_indel_file: ## OPEN file
        for line in SNV_indel_file: #Read in line by line
            line = line.rstrip("\n")
            columns = line.split("\t")
            if line[0] != "#": #Everything but header
                na12878 = columns[10] #Extracting individual NA12878
                if len(columns[3]) == 1 and len(columns[4]) == 1: #Specifying for the SNVs
                    snv_count += 1 #Add 1 to SNV count
                    if "0/0" in na12878: #If genotype is 0/0
                        SNV_homozygous_ref += 1 #Add to SNV homozygous reference count
                    elif "0/1" in na12878: #If genotype is 0/1
                        SNV_heterozygous += 1 #Add 1 to SNV het count
                    elif "1/1" in na12878: #If genotype is 1/1
                        SNV_homozygous_alternative += 1 #Add 1 to SNV homo alt count
                    elif "./." in na12878: #If genotype is missing
                        SNV_missing_genotypes += 1 #Add 1 to missing count
                        snv_count -= 1
                else: #Otherwise
                    indel_count += 1 #Add +1 to indel count
                    if "0/0" in na12878: #If 0/0
                        INDELS_homozygous_ref += 1 ##Add 1 to homo ref indel count
                    elif "0/1" in na12878: #if 0/1
                        INDELS_heterozygous += 1 #Add 1 to het count
                    elif "1/1" in na12878: #If 1/1
                        INDELS_homozygous_alternative += 1 #Add 1 to homo indel alt count
                    elif "./." in na12878: #Lastly, if ./. (missing)
                        INDELS_missing_genotypes += 1 #Add 1 to missing count
                        indel_count -= 1 #and subtract 1 from total indel count

## PRINT EVERYTHING

print("SNV Homozygous reference count: " + str(SNV_homozygous_ref))
print("SNV Heterozygous count: " + str(SNV_heterozygous))
print("SNV Homozygous alternative count: " + str(SNV_homozygous_alternative))
print("SNV Missing genotypes count: " + str(SNV_missing_genotypes))
print("Homozygous indel reference count: " + str(INDELS_homozygous_ref))
print("Heterozygous indel count: " + str(INDELS_heterozygous))
print("Homozygous alternative indel count: " + str(INDELS_homozygous_alternative))
print("Missing indel genotypes count: " + str(INDELS_missing_genotypes))
print("Indel count: " + str(indel_count))
print("SNV count: " + str(snv_count))

## COMBINE for total SNVs and Indel possibilities

Four_possibilities_SNVs = SNV_homozygous_ref+SNV_heterozygous+SNV_homozygous_alternative+SNV_missing_genotypes
Four_possibilities_Indels = INDELS_homozygous_ref+INDELS_heterozygous+INDELS_homozygous_alternative+INDELS_missing_genotypes

## Print 
print("Four possible SNV alleles total: " + str(Four_possibilities_SNVs))
print("Four possible indel alleles total: " + str(Four_possibilities_Indels))

