#!/usr/bin/env python3

"""make_seq.py prints a random sequence given a sequence length and nucleotide frequencies.  The random sequence will have the same nucleotide frequencies as the input nucleotide frequencies.

Usage: python3 make_seq.py <sequence_length> <a_freq> <c_freq> <g_freq> <t_freq>

<sequence_length> = Length of sequence to generate
<a_freq> = frequency of As
<c_freq> = frequency of Cs
<g_freq> = frequency of Gs
<t_freq> = frequency of Ts
"""

# Import modules 
import sys
import random

# sys.arg is a list containing 6 elements: the script name and 5 command line arguments
# Check that all 5 command line arguments were given. If not, print the documentation and exit.
if (len(sys.argv) != 6):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__) 

# Save the input arguments as variables
# By default, the command line arguments are saved as strings. Convert them to numeric types.
sequence_length = int(sys.argv[1])
a_freq = float(sys.argv[2])
c_freq = float(sys.argv[3])
g_freq = float(sys.argv[4])
t_freq = float(sys.argv[5])

# Check that frequencies add to 1. If not, exit the program 
if (abs(a_freq + t_freq + c_freq + g_freq - 1) > 1e-4):
	sys.exit("ERROR: Nucleotide frequencies do not add up to 1!")

## Part 4
### TODO Generate a random nucleotide sequence

#python3 make_seq.py <sequence_length> <A_freq> <C_freq> <G_freq>
# Initialize an empty string that nucleotides can be appended to
# Create a for loop that will be repeated <sequence_length> times
#     Generate a random decimal
#     Use if/else if/else logic to determine which nucleotide to add

Freq = [] #Initializes the list with all the nucleotides
for i in range(0,sequence_length): #loop for position 0 to len(sequence_length)
	choice = random.random() #Randomly draw a value 0-1
	if choice >= 1-a_freq: #if this decimal is between 1 and the designated value of a_freq in the python3 command, then append and print "A"
		Freq.append("A")
		print("A")
	elif choice >= 1-a_freq-t_freq: #otherwise, if the decimal is between 1 and the designated value of a_freq-t_freq in the python3 command, then append and print "T"
		Freq.append("T")
		print("T")
	elif choice >= 1-a_freq-t_freq-c_freq: #otherwise, if the decimal is between 1 and the designated value of a_freq-t_freq-c_freq in the python3 command, then append and print "T"
		Freq.append("C")
		print("C")
	else: #if the decimal falls outside of these above ranges, then print "G". 
		Freq.append("G")
		print("G")
