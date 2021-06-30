#!/usr/bin/env python3
"""
This script filters and normalizes RNA-Seq data, and produces several plots. 

Usage: python3 gene_expression.py raw_counts.txt
Arguments:
<list_of_samples>  = raw_counts.txt file (Req.)

"""

import sys
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import string
import operator

#TODO: Print out the doc string and exit if the number of input parameters is not correct

if len(sys.argv) != 2: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting two')

#######################
##Part 0 -- Functions##
#######################

# CPM counts per million
# Convert raw counts to counts per million (cpm)
# Raw count x[i,j] (from sample i, gene j)
# Total counts N[i] (from sample i)
# cpm[i,j] = (10^6)*x[i,j]/N[i]
def counts_per_million(dictionary, list_of_samples):
    # Initialize output dictionary cpm_dict
    cpm_dict = {}
    # Find N, the list of each sample's library size
    N = library_sizes(dictionary, list_of_samples)
    # Calculate cpm one gene at a time using list comprehension
    for k,v in dictionary.items():
        # Note: do NOT use 10^6, which equals 12
        # Use either 10**6 or 1000000 for 1 million
        # k is the key (gene name)
        # v is the value (raw counts of gene k)
        cpm_dict[k] = [(10**6)*x/n for x,n in zip(v,N)]
    return(cpm_dict)
# End of CPM function #

# Library size
# Calculate library size of each sample (e.g. sum of RNA-seq counts)
def library_sizes(dictionary, list_of_samples):
    # Get the total number of samples
    num_samples = len(list_of_samples)
    # Initialize N, a list to hold the total counts from each sample
    N = []
    # For loop to iterate over each sample 
    for i in range(num_samples):
        # Append a new float zero value for each sample (goes to index i)
        N.append(0.0)
        # For loop to iterate over each value in our dictionary
        for v in dictionary.values():
            # Get the count from the i index of this gene and add it to the total for sample i
            N[i] += v[i]
    # Return the list containing each sample's library size
    return(N)
# End of lbirary sizes function #

# TODO: Fill in the code for the translate_dictionary function.
# Translate dictionary
# Function to translate a dictionary from {gene:[list of counts by sample]} to {sample:[list of counts by gene]}
# The new dictionary will have one key for each sample and the value of each key will be a list of counts associated with that sample
# This function is used in another function called upper_quartile_norm()
# Each comment below should correspond to one line of code
def translate_dictionary(dictionary, list_of_samples):
    # Initialize empty dictionary called translated_dictionary
    translated_dictionary = {}
    #sort dictionary and extract each gene 
    each_gene = sorted(list(dictionary))
    #for loop to run through length of sample list
    for k in range(len(list_of_samples)):
        #extract sample name
        current_sample = list_of_samples[k]
        #new empty list for key with each sample name 
        translated_dictionary[current_sample] = []
        #for loop for every gene
        for i in each_gene:
            #identify count for ith gene in sample
            RNA_count = dictionary[i][k] 
            #append count for each sample in translated_dictionary
            translated_dictionary[current_sample].append(RNA_count)
    # Return the now translated dictionary
    return(translated_dictionary)

# TODO: Add comments to the code in the upper_quartile_norm function to explain what each line does
# Upper quartile normalization
# Compute the upper quartile normalization of raw counts
# Raw count x[i,j] (from sample i, gene j)
# D[i] value corresponding to 75th percentile of raw counts (from sample i)
# Mean of D, or the mean of all upper quartile means
# Upper quartile normalized count for sample i and gene j is (Mean of D)*x[i,j]/D[i]

def upper_quartile_norm(dictionary, list_of_samples):
    # initializing an empty dictionary
    upper_quartile_norm_dictionary = {}
    # define the number of items in the list_of_samples object, assign it to "num_samples"
    num_samples = len(list_of_samples)
    # Use the above translate dictionary function to translate, and assign it to "translated_dict_for_D", the now translated object
    translated_dict_for_D = translate_dictionary(dictionary, list_of_samples)
    # Initialize an empty list with straight brackets
    D = []
    # A for loop to iterate over each sample within num_samples
    for i in range(num_samples):
        # extract the current sample when running the loop, assign it to sample_name
        sample_name = list_of_samples[i]
        # Append 75th percentile of raw counts within the empty list 'D' for each sample
        D.append(np.percentile(translated_dict_for_D[sample_name],75))
    # Compute the arithmetic mean for all of D and assign to meanD
    meanD = np.mean(D)
    # Initialize for loop for two lists within the dictionary
    for k,v in dictionary.items():
        # for k, obtain the upper quartile normalized count for each respective sample/gene within the zipped lists
        upper_quartile_norm_dictionary[k] = [meanD*x/d for x,d in zip(v,D)]
    # return the dictionary with the upper quartile normalized counts
    return(upper_quartile_norm_dictionary)
	
# TODO: Write a function to calculate Fisher's Linear Discriminant (add comments, too!) for all genes in your count dictionary. Call your function fishers_linear_discriminant. Remember that functions should be able to work using different data sets, so make sure that your function would work with data from a different experiment with different numbers of before and after samples. 
# Return the top ten genes according to the highest FLD values
# Input to this function should be a count dictionary and two lists letting the function know which index values go with each group
# You will need to iterate over each gene in the dictionary and calculate the FLD of each gene
# FLD of a gene = ((m1-m2)^2)/((s1)^2 + (s2)^2)
# m1 = mean of the first group
# m2 = mean of the second group
# s1 = standard deviation of the first group
# s2 = standard deviation of the second group

def fishers_linear_discriminant(dictionary, list_one, list_two):
	#initialize FLD dictionary for each gene
	FLD_dictionary = {}
	#for loop for keys and values in dictionary 
	for k,v in dictionary.items():
	#empty list for group 1
		GroupOne = []
	#empty list for group 2
		GroupTwo = []
	#Initiate for loop for each sample in list_one
		for x in list_one:
	#append value into  empty list GroupOne
			GroupOne.append(v[x])
	#Initiate for loop for each sample in list_two
		for y in list_two:
	#append value into empty list GroupTwo
			GroupTwo.append(v[y])
		#calculate mean and SD
		m1 = np.mean(GroupOne)
		s1 = np.std(GroupOne)
		m2 = np.mean(GroupTwo)
		s2 = np.std(GroupTwo)
		#Calculate FLD and store next to gene
		current_FLD = ((m1-m2)**2)/((s1)**2 + (s2)**2)
		#Return the dictionary with each FLD, put keys/values in!
		FLD_dictionary[k] = current_FLD
	return(FLD_dictionary)
  
		
####################
##End of functions##
####################

# These first lines of code will get the data imported and in the right format for the rest of the homework
# You need to do the rest of the work starting from Part 1 -- Data filtering

#open the data file
data_file = open(sys.argv[1],"r")
#first line of data file contains sample names -- store in a list
sample_list = data_file.readline().strip().split()[1:]
#initialize raw count dictionary
raw_counts_dict = {}
#add each gene and expression values to dictionary
for line in data_file:
    line_list = line.strip().split() #split line into list at whitespace, strip removes leading and trailing whitespace
    gene = line_list[0] #name of gene is the first thing in the list
    expression_values = [int(float(v)) for v in line_list[1:]] #values of gene expression follow the gene name
    raw_counts_dict[gene] = expression_values #add keys and values to the dictionary {gene:expression}
#close the data file
data_file.close()

############################
##Part 1 -- Data filtering##
############################


# Filter out genes with zero expression in all samples

#print(len(raw_counts_dict)) ##full dict has 27012 genes

#Extracting a dictionary using dictionary comprehension to include only genes that have values that all DO NOT sum to zero.
#Source: https://stackoverflow.com/questions/49938985/filter-out-entries-in-dictionary-based-on-the-sum-of-the-keys
NoZeros = {k: v for k, v in raw_counts_dict.items() if sum(v) != 0}  

# Printing how many genes are in the filtered dictionary
#use str() to combine with a more informative string about the value itself
print("Number of genes without all zero expression counts: " + str(len(NoZeros)))

# Filter out genes with 20 or more samples with cpm < 1

#Convert raw counts to CPM and assign to "OneCPM" 
OneCPM = counts_per_million(NoZeros,sample_list)
#Initiate several empty dictionaries:
Good_Genes_FilterTwo = {} #Will hold genes that pass second filter
Bads_Genes_FilterTwo = {} #Will hold genes that do not pass second filter
Good_Genes_FilterOne = {} #Will hold genes that passed filter one
#Initial counts for second filter dictionaries will be 0
Good_Genes_FilterTwo_Count = 0
Bad_Genes_FilterTwo_Count = 0
# Assigning cpm_dict to counts_per_million of OneCPM
OneCPM = counts_per_million(NoZeros, sample_list)
#for loop for keys and values within NoZeros using .items() for tuples
for k,v in NoZeros.items():
    #The initial value for the failed CPM list should be zero
    CPM_Failed = 0
    #for loop to iterate through OneCPM
    for CPM_Value in OneCPM[k]:
        #if the value of a CPM is less than one:
        if CPM_Value < 1:
            #Then increase the count by one
            CPM_Failed += 1
    #If the value of a CPM is greater or equal than 20
    if CPM_Failed >= 20:
        #Assign into the bad genes filter
        Bads_Genes_FilterTwo[k] = v
        #Add +1 to the counter to update 
        Bad_Genes_FilterTwo_Count += 1
        #Else
    else:
        #Otherwise, add to good genes filter
        Good_Genes_FilterTwo[k] = v
        #Add +1 to the counter for good genes
        Good_Genes_FilterTwo_Count += 1
#Print number of genes that make it past this filter with a combined string
print("Number of filtered genes excluding those with 20 or more samples with CPM < 1: " + str(Good_Genes_FilterTwo_Count))
 
################################
##Part 2 -- Data visualization##
################################
#Create a tuple list of all genes based on the summed values for each sample- NOTE: these are without sample identifiers, thus I am assuming the order is kept the same-- unsure if this will be OK in this case.
summed = tuple(sum(x) for x in zip(*Good_Genes_FilterTwo.values()))

#Turn these summed values in a dictionary where each sample has a unique numeric identifier
#https://stackoverflow.com/questions/47405022/python-create-dictionary-with-list-index-number-as-key-and-list-element-as-valu
#Initialize empty dictionary to hold tupple values with unique numeric identifier
Unique_Summed={}
#For loop for key and value to enumerate() each observation within the tupple list and assign a value 
for index,value in enumerate(summed):
    Unique_Summed[index+1]=value #add one for each value for 1-40

#Now we need to replace numbers with sample names. We can do this by matching the dictionary to the order of the sample_list
SampleNamesSampleCounts=dict(zip(sample_list,list(Unique_Summed.values()))) 

#Now we shall plot. 
plt.figure(figsize=(10, 10)) #Define plot dimensions
plt.bar(list(SampleNamesSampleCounts.keys()), SampleNamesSampleCounts.values(), color='lightskyblue') #In sky blue, make a barplot for the value in each key
plt.ylabel('Total Read Count (million)') #Define label for y-axis
plt.xlabel('Sample Number (Identifier)') #Define label for x-axis
plt.title('Total RNA-Seq reads counts after filtering') #Define main title
plt.xticks(fontsize=9, rotation=90) #Make font a bit smaller, and rotate the labels 90 degrees to fit many samples.
plt.xticks(np.arange(0, 40, 1)) #Arrange the x-tick labels appropriately (40 ticks, one label per tick spaced across the plot)
plt.savefig('library_size.png') #Save the plot

#Get range
#Source: https://stackoverflow.com/questions/18595686/how-do-operator-itemgetter-and-sort-work
#Identify the maximum value sample ('max')
print("Sample with highest read count after filtering: " + str(max(SampleNamesSampleCounts.items(), key=operator.itemgetter(1))[0]))
#Max is sample After_1 = 17301749
#Identify the minimum value sample ('min')
print("Sample with lowest read count after filtering: " + str(min(SampleNamesSampleCounts.items(), key=operator.itemgetter(1))[0]))
#Min is sample After_19 = 11464513

################################
##Part 3 -- Data normalization##
################################

# Normalize count data left after filtering steps
#run the upper_quartile_norm function on the genes that passed the second filter
upper_quartile_norm_results = upper_quartile_norm(Good_Genes_FilterTwo, sample_list)
#create summed count value tupple subsequent to normalization
summed_normalized = tuple(sum(x) for x in zip(*upper_quartile_norm_results.values()))
#Turn these summed values in a dictionary where each sample has a unique numeric identifier-- similar to second filter protocol
#Source: https://stackoverflow.com/questions/47405022/python-create-dictionary-with-list-index-number-as-key-and-list-element-as-valu
#Initialize empty dictionary
dict_normalized={}
#Quick for loop to enumerate each value in dictionary uniquely
for index,value in enumerate(summed_normalized):
    dict_normalized[index+1]=value #Add one for 1-40

#Replace numbers with sample names 
NormalizedSampleCountsNames=dict(zip(sample_list,list(dict_normalized.values()))) 

#Plot figure
plt.figure(figsize=(10, 10)) #Define plot dimensions
plt.bar(list(NormalizedSampleCountsNames.keys()), NormalizedSampleCountsNames.values(), color='coral') #In coral, make a barplot for the value in each key
plt.ylabel('Total Read Count (million)') #Define label for y-axis
plt.xlabel('Sample Number (Identifier)') #Define label for x-axis
plt.title('Total RNA-Seq reads counts after filtered normalization') #Main title
plt.xticks(fontsize=9, rotation=90) #Font size and tick label rotation
plt.xticks(np.arange(0, 40, 1)) #Arrange the x-tick labels appropriately (40 ticks, one label per tick spaced across the plot)
plt.savefig('library_size_normalized.png') #Save the plot

##Get range
#Print max normalized count sample with combined string
#Max is sample After_13
print("Sample with highest read count after filtered normalization: " + str(max(NormalizedSampleCountsNames.items(), key=operator.itemgetter(1))[0]))
#Print min normalized count sample with combined string
#Min is sample After_19
print("Sample with lowest read count after filtered normalization: " + str(min(NormalizedSampleCountsNames.items(), key=operator.itemgetter(1))[0]))
##############################
##Part 4 -- Data exploration##
##############################

# Calculate Fisher's Linear Discriminant for each gene basd on the normalized count data
#Calculate FLD for each gene, use range() to select between Before and After groups
FishersValues = fishers_linear_discriminant(upper_quartile_norm_results, range(20), range(20,40))
#Print to make sure everything looks OK
#print(FishersValues)
#Sort the list from highest to lowest, selecting top 10 genes
#Source: https://stackoverflow.com/questions/7197315/5-maximum-values-in-a-python-dictionary
FishersValuesTopTen = sorted(FishersValues, key=FishersValues.get, reverse=True)[:10]

#Print the top ten genes with corresponding FLD scores
#print(FishersValuesTopTen)
print("Top 10 genes based on FLD values:")
#for loop to match up these top ten gene names with actual FLD scores
for i in FishersValuesTopTen:
	print(i + '=' + str(FishersValues[i]))
	
## Pick a gene to explore further and plot the mean expression level for the before and after groups (save file as mean_expression.png)
#MYLK4
#https://problemsolvingwithpython.com/06-Plotting-with-Matplotlib/06.07-Error-Bars/
MYLK4 = ['MYLK4'] #Select MYLK4
#Quick for loop to pull MYLK4 data in the normalized count dataset
for normalized_count in MYLK4:
	MYLK = upper_quartile_norm_results[normalized_count] #Assign to MYLK
#Before and After groups as numerical data
MYLK_Before = np.array(MYLK[0:20]) 
MYLK_After = np.array(MYLK[20:40])

#Calculate means
MYLK_Before_mean = np.mean(MYLK_Before)
MYLK_After_mean = np.mean(MYLK_After)

#Calculate SDs
MYLK_Before_SD = np.std(MYLK_Before)
MYLK_After_SD = np.std(MYLK_After)

labels = ['MYLK4_Before', 'MYLK4_After'] #Define x-axis tick labels
x_pos = np.arange(len(labels)) #Define space between ticks
CTEs = [MYLK_Before_mean, MYLK_After_mean] #Combine means together for plotting
error = [MYLK_Before_SD, MYLK_After_SD] #Combine SDs together for plotting

fig, ax = plt.subplots() #Define plot
ax.bar(x_pos, CTEs, #Define the error bars with attributes (e.g., black, centered, transparency)
       yerr=error,
       align='center',
       alpha=0.5,
       ecolor='black',
       capsize=10)
ax.set_ylabel('Average Read Count') #Provide y-axis label
ax.set_xticks(x_pos) #Assign x-axis tick positions
ax.set_xticklabels(labels) #Controls appearance of tick labels
ax.set_title('MYLK4 Before and After Group Counts (normalized)') #Main title
ax.yaxis.grid(True) #Gridded

# Save the figure and show
plt.tight_layout() #For export layout
plt.savefig('mean_expression.png') #Export plot

### END SCRIPT ###