#!/usr/bin/env python3

""" 
nuc_count.py counts nucleotides in a fasta file

Usage: python3 nuc_count.py <fasta>

<fasta> = path to a fasta file
""" 

# Import modules
import sys
import re
# sys.arg is a list containing 2 elements: the script name and 1 command line argument
# Check that all the command line argument was given. If not, print the documentation and exit.
if (len(sys.argv) != 2):
    sys.exit(__doc__) 

# Save the input arguments as variables
fasta = sys.argv[1]

# Initialize a nucleotide string
nucleotides = ""

# Load the fasta sequence
# NOTE: this script assumes there is only *one* sequence in the fasta file
# Open the fasta file
with open(fasta) as f:
    # For each line in the file
    for line in f:
        # Skip lines starting with ">"
        if not line.startswith(">"):
            # Add each line to the nucleotide string
            nucleotides += line.rstrip()

# Make the nucleotide string all capital letters 
nucleotides = nucleotides.upper()

# Count the nucleotides and print output
num_a = nucleotides.count('A')
num_c = nucleotides.count('C')
num_g = nucleotides.count('G')
num_t = nucleotides.count('T')
#num_n = nucleotides.count('N')

print ("Nucleotide Raw Counts")
print ("A: ", num_a)
print ("C: ", num_c)
print ("G: ", num_g)
print ("T: ", num_t)
#print ("N: ", num_n)


## Part 3
### TODO Print out the frequencies for each nucleotide in alphabetical order
base_sum = num_a + num_c + num_g + num_t #sum the total As, Ts, Cs, and Gs and assign to base_sum
print ("Base Sum", base_sum) #Print base sum
print ("Frequency of A:", round(num_a/base_sum, 2)) #Print frequency of each nucleotide and round
print ("Frequency of C:", round(num_c/base_sum, 2))
print ("Frequency of G:", round(num_g/base_sum, 2))
print ("Frequency of T:", round(num_t/base_sum, 2))

# Count the nucleotides and print output
num_AA = nucleotides.count('AA')
num_AG = nucleotides.count('AG')
num_CA = nucleotides.count('CA')
num_CG = nucleotides.count('CG')
num_GA = nucleotides.count('GA')
num_GG = nucleotides.count('GG')
num_TA = nucleotides.count('TA')
num_TG = nucleotides.count('TG')
num_AC = nucleotides.count('AC')
num_AT = nucleotides.count('AT')
num_CC = nucleotides.count('CC')
num_CT = nucleotides.count('CT')
num_GC = nucleotides.count('GC')
num_GT = nucleotides.count('GT')
num_TC = nucleotides.count('TC')
num_TT = nucleotides.count('TT')

# Count the number of dinucleotides -- this will print the raw number
print ("Dinucleotides Raw Counts")
print ("AA: ", num_AA)
print ("AG: ", num_AG)
print ("CA: ", num_CA)
print ("CG: ", num_CG)
print ("GA: ", num_GA)
print ("GG: ", num_GG)
print ("TA: ", num_TA)
print ("TG: ", num_TG)
print ("AC: ", num_AC)
print ("AT: ", num_AT)
print ("CC: ", num_CC)
print ("CT: ", num_CT)
print ("GC: ", num_GC)
print ("GT: ", num_GT)
print ("TC: ", num_TC)
print ("TT: ", num_TT)

## Part 3
# Print the frequencies for each dinucleotide pair
#Sum total dinucleotides below
di_base_sum = num_AA + num_AG + num_CA + num_CG + num_GA + num_GG + num_TA + num_TG + num_AC + num_AT + num_CC + num_CT + num_GC + num_GT + num_TC + num_TT
print ("Dinucleotide Base Sum", di_base_sum) #Print the dinucleotide base sum
print ("Frequency of AA:", round(num_AA/di_base_sum, 2)) #Print each dinucleotide as a frequency (dinucleotide/dinucelotide total) and round
print ("Frequency of AG:", round(num_AG/di_base_sum, 2))
print ("Frequency of CA:", round(num_CA/di_base_sum, 2))
print ("Frequency of CG:", round(num_CG/di_base_sum, 2))
print ("Frequency of GA:", round(num_GA/di_base_sum, 2))
print ("Frequency of GG:", round(num_GG/di_base_sum, 2))
print ("Frequency of TA:", round(num_TA/di_base_sum, 2))
print ("Frequency of TG:", round(num_TG/di_base_sum, 2))
print ("Frequency of AC:", round(num_AC/di_base_sum, 2))
print ("Frequency of AT:", round(num_AT/di_base_sum, 2))
print ("Frequency of CC:", round(num_CC/di_base_sum, 2))
print ("Frequency of CT:", round(num_CT/di_base_sum, 2))
print ("Frequency of GC:", round(num_GC/di_base_sum, 2))
print ("Frequency of GT:", round(num_GT/di_base_sum, 2))
print ("Frequency of TC:", round(num_TC/di_base_sum, 2))
print ("Frequency of TT:", round(num_TT/di_base_sum, 2))

## Part 5
### TODO Use overlapping windows to count the dinucleotides. See the assignment for more information on overlapping windows.
### TODO Print the frequencies for each of dinucleotides in alphabetical order

### This stackexchange example was very helpful https://bioinformatics.stackexchange.com/questions/2043/how-can-i-do-an-overlapping-sequence-count-in-biopython 

AA_overlap = len(re.findall(r'(?=(AA))',nucleotides)) #Calculate overlapping window nucleotide lengths using "re.findall" from re module
AG_overlap = len(re.findall(r'(?=(AG))',nucleotides))
CA_overlap = len(re.findall(r'(?=(CA))',nucleotides))
CG_overlap = len(re.findall(r'(?=(CG))',nucleotides))
GA_overlap = len(re.findall(r'(?=(GA))',nucleotides))
GG_overlap = len(re.findall(r'(?=(GG))',nucleotides))
TA_overlap = len(re.findall(r'(?=(TA))',nucleotides))
TG_overlap = len(re.findall(r'(?=(TG))',nucleotides))
AC_overlap = len(re.findall(r'(?=(AC))',nucleotides))
AT_overlap = len(re.findall(r'(?=(AT))',nucleotides))
CC_overlap = len(re.findall(r'(?=(CC))',nucleotides))
CT_overlap = len(re.findall(r'(?=(CT))',nucleotides))
GC_overlap = len(re.findall(r'(?=(GC))',nucleotides))
GT_overlap = len(re.findall(r'(?=(GT))',nucleotides))
TC_overlap = len(re.findall(r'(?=(TC))',nucleotides))
TT_overlap = len(re.findall(r'(?=(TT))',nucleotides))

print ("Dinucleotides Overlapping Raw Counts") #Print number of overlapping dinucleotides for each pair
print("AA Overlap: ", AA_overlap)
print("AG Overlap: ", AG_overlap)
print("CA Overlap: ", CA_overlap)
print("CG Overlap: ", CG_overlap)
print("GA Overlap: ", GA_overlap)
print("GG Overlap: ", GG_overlap)
print("TA Overlap: ", TA_overlap)
print("TG Overlap: ", TG_overlap)
print("AC Overlap: ", AC_overlap)
print("AT Overlap: ", AT_overlap)
print("CC Overlap: ", CC_overlap)
print("CT Overlap: ", CT_overlap)
print("GC Overlap: ", GC_overlap)
print("GT Overlap: ", GT_overlap)
print("TC Overlap: ", TC_overlap)
print("TT Overlap: ", TT_overlap)

# Sum the total dinucleotides
Di_Base_Sum_Overlapping = AA_overlap + AG_overlap + CA_overlap+ CG_overlap +GA_overlap +GG_overlap +TA_overlap +TG_overlap +AC_overlap +AT_overlap +CC_overlap +CT_overlap +GC_overlap +GT_overlap +TC_overlap +TT_overlap 
# Print each dinucleotide combination as a frequency (e.g., overlapping AA/total, round to two decimal places)
print ("Dinucleotide Overlapping Frequencies") #Print overlapping dinucleotide frequenciesfor each pair
print ("AA Overlapping Frequency:", round(AA_overlap/Di_Base_Sum_Overlapping, 2))
print ("AG Overlapping Frequency:", round(AG_overlap/Di_Base_Sum_Overlapping, 2))
print ("CA Overlapping Frequency:", round(CA_overlap/Di_Base_Sum_Overlapping, 2))
print ("CG Overlapping Frequency:", round(CG_overlap/Di_Base_Sum_Overlapping, 2))
print ("GA Overlapping Frequency:", round(GA_overlap/Di_Base_Sum_Overlapping, 2))
print ("GG Overlapping Frequency:", round(GG_overlap/Di_Base_Sum_Overlapping, 2))
print ("TA Overlapping Frequency:", round(TA_overlap/Di_Base_Sum_Overlapping, 2))
print ("TG Overlapping Frequency:", round(TG_overlap/Di_Base_Sum_Overlapping, 2))
print ("AC Overlapping Frequency:", round(AC_overlap/Di_Base_Sum_Overlapping, 2))
print ("AT Overlapping Frequency:", round(AT_overlap/Di_Base_Sum_Overlapping, 2))
print ("CC Overlapping Frequency:", round(CC_overlap/Di_Base_Sum_Overlapping, 2))
print ("CT Overlapping Frequency:", round(CT_overlap/Di_Base_Sum_Overlapping, 2))
print ("GC Overlapping Frequency:", round(GC_overlap/Di_Base_Sum_Overlapping, 2))
print ("GT Overlapping Frequency:", round(GT_overlap/Di_Base_Sum_Overlapping, 2))
print ("TC Overlapping Frequency:", round(TC_overlap/Di_Base_Sum_Overlapping, 2))
print ("TT Overlapping Frequency:", round(TT_overlap/Di_Base_Sum_Overlapping, 2))


