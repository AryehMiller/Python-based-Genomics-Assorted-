#!/usr/bin/env python3
"""Scan a DNA sequence to find putative binding sites

Usage: python3 scan_sequence.py <scoring_matrix> <sequence_file> <threshold>

Args:
    scoring_matrix = Path to scoring matrix. The rows of the matrix correspond 
       to A, C, G, and T, and the columns correspond to positions
    sequence_file = Path to DNA sequence file.
    threshold = numeric value (integer) that sets threshold to filter report hits by
"""
import sys

# Helper dictionary for reverse complementation
reverse_comp = { 'A':'T', 'C':'G', 'G':'C', 'T':'A' }

###############################################################
# Begin functions
###############################################################

def create_scoring_matrix_from_file(matrix_file):
    """
    for reading and returning a scoring matrix by pulling out the columns for downstream usage
    Args: matrix_file -- a file with columns corresponding to A,C,G, and T
    """
    file_data = [ x.split() for x in open(matrix_file).readlines() ]
    scoring_matrix = [ dict(A=float(a_score), C=float(c_score), G=float(g_score), 
        T=float(t_score)) for a_score, c_score, g_score, t_score in 
        zip(file_data[0], file_data[1], file_data[2], file_data[3])
    ]

    return scoring_matrix


def read_sequence_from_file(file_name):
    """
    this function reads in the sequence from the sequence file and strips the lines
    Args: file_name -- a DNA sequence file with a string of bases
    """
    return open(file_name).readlines()[0].strip()



def score_with_matrix(subseq, matrix):
    """
    this function returns the sum of scores for a given sequence of base pairs using the scoring matrix
    Args: subseq -- a partial sequence, matrix -- scoring matrix with columns corresponding to A,C,G, and T
    """
    return sum([ score[ base ] for score, base in zip(matrix, subseq)])


def get_reverse_complement(in_string):
    """
    This function takes a nucleotide sequence and returns the reverse compliment
    Args: in_string -- a sequence of nucleotides
    """
    return (''.join([ reverse_comp[x] for x in in_string[::-1] ]))

###############################################################
# End functions
###############################################################

###############################################################
# Begin main script
###############################################################

# Check the correct number of command line arguments
if(len(sys.argv)!= 4):
    sys.exit(__doc__)
    
#assign files
score_matrix_file = sys.argv[1]
sequence_file = sys.argv[2]
threshold = float(sys.argv[3])

# TODO: explain what this code does
score_matrix = create_scoring_matrix_from_file(score_matrix_file) #running the function 'create_scoring_matrix_from_file' on the first argument 'score_matrix_file' to have a readable and workable scoring matrix file
motif_width = len(score_matrix) #get the number of items in score_matrix and assign number to motif_width
search_sequence = read_sequence_from_file(sequence_file) #read in the sequence file and assign to search_sequence
#Generate the reverse complement
reverse_sequence = get_reverse_complement(search_sequence) #take the reverse complement

# Calculate the number of matrix 'windows' for calculating sequence scores
last_index = len(search_sequence) - motif_width + 1

#Creates a list that joins orientation, positions, and scores for all occurrences
forward_hit_list = [ (i, score_with_matrix(search_sequence[i:i+motif_width], score_matrix)) for i in range(last_index) ] #forward
reverse_hit_list = [ (i, score_with_matrix(reverse_sequence[i:i+motif_width], score_matrix)) for i in range(last_index) ] #reverse

#for loop for forward and reverse to run through the list and find observations that meet the threshold (greater or equal to "threshold")
#Source: https://www.geeksforgeeks.org/python-find-the-tuples-containing-the-given-element-from-a-list-of-tuples/
forward_hit_list = [item for item in forward_hit_list 
          if item[1] >= threshold] #filter for greater than or equal to the specified argumental threshold
          
reverse_hit_list = [item for item in reverse_hit_list 
          if item[1] >= threshold]  #filter for greater than or equal to the specified argumental threshold
              
# Forward
if len(forward_hit_list) == 0: #Checking the length of the forward and reverse sequences, specifically if they are equal to 0
    print("No threshold-exceeding hits found in the forward direction!") #if true (no observations meet condition of exceeding threshold), then print that
else: #otherwise, continue on
    print("orientation\tsequence\tposition\tscore") #print column headers
    for hit in forward_hit_list: #for loop for each position in 
        print("forward\t{sequence:s}\t{position:d}\t{score:.2f}".format(position=hit[0], sequence=search_sequence[hit[0]:hit[0]+motif_width],score=hit[1])) #print forward entries table-- include the sequence=search_sequence to loop in the actual corresponding sequence
        
# Reverse complement
if len(reverse_hit_list) == 0: #Checking the length of the forward and reverse sequences, specifically if they are equal to 0
    print("No threshold-exceeding hits found in reverse direction!") #if true, then print that
else: #otherwise, continue on
    print("orientation\tsequence\tposition\tscore") #print column headers
    for hit in reverse_hit_list: #for loop for each position
        print("reverse\t{sequence:s}\t{position:d}\t{score:.2f}".format(position=hit[0], sequence=reverse_sequence[hit[0]:hit[0]+motif_width],score=hit[1])) #print reverse entries table-- include the sequence=search_sequence to loop in the actual corresponding sequence


