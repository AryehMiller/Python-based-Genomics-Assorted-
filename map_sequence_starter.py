#!/usr/bin/env python3
"""
Align sequencing reads to a set of genes.

Usage: map_sequence_starter.py <cDNA_fasta> <seq_reads_file>

Arguments:
    <cDNA_fasta>     = Path to cDNA fasta file.  (Required.)
    <seq_reads_file> = Path to sequencing reads file.  (Required.)
"""

import sys
#Source for the below code: https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
def reverse_complement(seq_reads_file):
    #This function takes the reverse complement of a sequence https://gist.github.com/hurshd0/3d54301cdb052cc4c103a2d118a9c7b9
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N':'N'} #RC nucleotides
    reverse_c_string = seq_reads_file [::-1] #reverse the string
    reverse_comp = [] #make empty list for the reverse complemented sequence
    seq_list = list(seq_reads_file) #create a list out of the seq_reads_file
    for base in reverse_c_string: #initiate for loop for each base in the RC string
        #append each complementary base from the dictionary on line 16 in the empty list
        reverse_comp.append(complement_dict[base])
    #return the reverse complement string and join RC bases together
    return ''.join(reverse_comp)

def read_cDNA_file_to_dict(filename):   
    #initialize the dictionary
    cDNA_dict = {}
    #open the cDNA file
    with open(filename) as file_one:
        #loop through each line in file
        for line in file_one:
            #remove newline using .rstrip()
            line = line.rstrip()
            # get gene name
            #since this is a fasta file, if line begins with ">" denotes a new cDNA sequence (gene)
            if line [0] == '>':
                # gene names are within "|", so use .split to break up words in string
                gene_list = line.split("|")
                #select gene name
                gene_name = gene_list[1] 
            #if the first character in a line is not >, then it is the sequence
            else:
                #read in the sequence in uppercase to make sure there aren't any issues further downstream when it comes to matching nucleotides
                line = line.upper()
            #add name and sequence in dictionary
            cDNA_dict[gene_name] = line
    #return dictionary    
    return cDNA_dict 
    
def create_twenty_five_mer_dict(cDNA_dict):
    #This function creates a dictionary of 25mer sequences from the cDNA sequences.  The dictionary keys are the 25mer sequences and the values are the corresponding gene names.   
    #Source for some of the below code: https://www.biostars.org/p/268007/
    #initialize dictionary
    twenty_five_mer_dict = {}
    #loop through cDNA dictionary
    for key in cDNA_dict:
        #get sequence that corresponds to the key
        sequence = cDNA_dict[key]
        #move through sequence, grabbing 25 mers and place in dictionary
        #note: if you want to change the k-mer length, simply change 25 to the desired length
        for i in range(0, len(sequence)-25+1):
            #the object twenty_five_mer is the sequence from i to i+25
            twenty_five_mer = sequence[i:i+25]
            #add 25mer and gene name to twenty_five_mer_dict dictionary
            twenty_five_mer_dict[twenty_five_mer] = key
    #return dictionary
    return twenty_five_mer_dict
    
##############main loop###################

#parse command line
#check that the correct number of arguments were given
if (len(sys.argv) != 3):
    sys.exit(__doc__)

cDNA_file = sys.argv[1]
seq_reads_file = sys.argv[2]

#read in cDNAs
cDNA_dict = read_cDNA_file_to_dict(cDNA_file)
print("Read in the cDNAs")


twenty_five_mer_dict = create_twenty_five_mer_dict(cDNA_dict)
print("Created dictionary")


count_dict = {}
for key in cDNA_dict:
    count_dict[key] = 0


fh = open(seq_reads_file, "r")

for line in fh:

    line = line.rstrip()
    line = line.upper() 
    if line in twenty_five_mer_dict: 
        if twenty_five_mer_dict[line] in count_dict:
            count_dict[twenty_five_mer_dict[line]] = count_dict[twenty_five_mer_dict[line]] + 1

    else: 
        rc_seq = reverse_complement(line)
        if rc_seq in twenty_five_mer_dict: 
            if twenty_five_mer_dict[rc_seq] in count_dict:
                count_dict[twenty_five_mer_dict[rc_seq]] = count_dict[twenty_five_mer_dict[rc_seq]] + 1

fh.close()

print("Name\tReads\tReads per BP")
for name, count in sorted(count_dict.items()):
    print(name + "\t" + str(count) + "\t" + str(float(count)/len(cDNA_dict[name])))
