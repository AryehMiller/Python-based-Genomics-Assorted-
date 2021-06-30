#!/usr/bin/env python3
"""Wright fisher model simulation

Usage: python3 wrightfisher.py <population size> <max # generations> <fitness value> <model>
python3 wrightfisher.py 100 1000 3 dominant
Args: <population size> an integer specifying the number of individuals in the simulated population
	  <max # generations > the number of generations in the simulation
	  < fitness value > the specified fitness value which is simulated
	  < model > either a "dominant" or "recessive" mode of genetic inheritance to model

"""
import sys
import csv
import re
import io
import os
import matplotlib.pyplot as plt
import numpy as np
import random

## Quit if input parameters aren't correct
if len(sys.argv) != 5: #if the arguments isn't correct (2), then print the doc string and the error indicating incorrect input parameters
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting four input parameters-- models is "dominant" or "recessive"')

## LOAD IN MODEL PARAMETERS

#Number of individuals in population
population_size = int(sys.argv[1]) ## Need Int
#Generations
Generations = int(sys.argv[2]) #Need Int
#Mutant allele fitness value (Bio5488)
fitness = float(sys.argv[3]) ## Need FLOAT since we're importing a numeric value for downstream use
#Mode of inheritance (either "dominant" or "recessive")
model = sys.argv[4]

## MODEL

## MODEL STRUCTURE
# 1) Simulate a biallelic population in a generation, a mere single mutant individual is present (either recessive or dominant, though)
# 2) Calculate frequency of alleles in that population, if the frequency of an allele (homozygous reference, homozygous mutant, heterozygous) is fixed, exit simulation and print fixation count
# 3) Repeat for number of simulated times
# 4) Calculate the proportion of time fixed


## We shall start our count of iterations until fixation
prop_fixated = 0
## First specify the three genotypes we might expect-- homozygous reference, heterozygous, and homozygous alternate (mutant)
quant_genotypes = ["00", "01", "11"]
## Placeholder
FitBio5488 = fitness
## Initiate for loop to run the simulation 100 times
for i in range(100):
	#Initiate population of homozygotes respective to the population size
    genotypes = [0] * population_size
    #Initiate AllFitnesses of 1 for all of the above individuals (homozygous references)
    AllFitnesses = [1] * population_size 
	#Use the random module to select a random individual
    Bio5488 = random.randrange(population_size)
    #Mutation occurs, now heterozygous individual
    genotypes[Bio5488] = 1
    #Here, I am stating that if the input parameter is "dominant" (mode of inheritance), then the individual will be heterozygous
    if model == "dominant":
    	#Mutation
        AllFitnesses[Bio5488] = FitBio5488
    ## MAKE SEVERAL OPEN LISTS
    #Reference allele frequency
    R_frequency = []
    #Homozygous reference
    RR = []
    #Heterozygous
    GR = []
    #Homozygous mutant
    GG = []
    R_frequency.append(sum(genotypes)/(2*population_size)) #Append reference allele freq and sum and divide them by 2 times the population size for prob distrib
    RR.append(genotypes.count(0)/population_size) #Append RR freq
    GR.append(genotypes.count(1)/population_size) #Append GR freq
    GG.append(genotypes.count(2)/population_size) #Append GG freq
    #Initiate for loop to loop through generations
    for i in range(Generations):
    	#Have AllFitnesses sum to 1
        FitProbabilities = [fitness/sum(AllFitnesses) for fitness in AllFitnesses]
        #Initiate list for the novel genotypes
        NovelGenotypes = []
        #Initiate list for novel fitnesses
        NovelFitnesses = []
        #Initiate for loop for subsequent generations
        for j in range(population_size):
            #Pick the first allele
            #Pick the genotype of a random individual
            #Source for weighted prob: https://www.geeksforgeeks.org/how-to-get-weighted-random-choice-in-python/
            #Randonly choose first genotypes/allele using a weighted probability format, assign to Parent1 (parent 1)
            Parent1 = genotypes[np.random.choice(population_size, p = FitProbabilities)]
            #Get allele
            Parent1Allele = random.randrange(2)
            #Binary allele (integer)
            Parent1Allele = int(quant_genotypes[Parent1][Parent1Allele])
			#SELECT SECOND ALLELE FROM OTHER PARENT
			#Pick the second individual
			#Pick the second allele
            #Randonly choose second genotypes/allele using a weighted probability format, assign to Parent2 (parent 2)
            Parent2 = genotypes[np.random.choice(population_size, p = FitProbabilities)]
            #Pick an allele from that genotyp
            Parent2Allele = random.randrange(2)
            #Binary allele (integer)
            Parent2Allele = int(quant_genotypes[Parent2][Parent2Allele])
            #Combine to ascertain the genotype of the subsequent generation
            NovelIndividualGenotype = Parent1Allele + Parent2Allele
            #Append it
            NovelGenotypes.append(NovelIndividualGenotype)
            #If statement-- if model is dominant and the genotype is 1, or the genotype is 2, then assign to the mutant fitness
            if (model == "dominant" and NovelIndividualGenotype == 1) or NovelIndividualGenotype == 2:
                NovelFitnesses.append(FitBio5488)
            #Otherwise, append 1
            else:
                NovelFitnesses.append(1)
        #Equate
        #Same as above block at beginning
        genotypes = NovelGenotypes
        AllFitnesses = NovelFitnesses
        R_frequency.append(sum(genotypes)/(2*population_size))
        RR.append(genotypes.count(0)/population_size)
        GR.append(genotypes.count(1)/population_size)
        GG.append(genotypes.count(2)/population_size)
		## NOW we can check for fixation
		#If the genotype has fixated throughout pop, break
        if genotypes == [2] * population_size:
            prop_fixated += 1
            break
        elif genotypes == [0] * population_size:
            break
print("Generations fixed: ", prop_fixated)

## PLOTTING

## Select last gen and plot
G_frequency = [1 - i for i in R_frequency]
## Use different line style (--) to plot the R and G frequencies
plt.plot(R_frequency,range(len(G_frequency)), linestyle='--', color='lightskyblue', label='R frequency')
plt.plot(G_frequency,range(len(G_frequency)), linestyle='--', color='red', label='G frequency')
#Genotype frequencies
plt.plot(GG, range(len(R_frequency)), color = "black", label = "GG frequency")
plt.plot(RR, range(len(R_frequency)), color = "blue", label = "RR frequency")
plt.plot(GR, range(len(R_frequency)), color = "green", label = "GR frequency")
#Plot labels
plt.xlabel('Generation slice')
plt.ylabel('Genotype frequencies')
#If statement for labeling based on model
if model == "dominant":
	plt.title('Simulated Wright-Fisher Model: Dominant')
elif model == "recessive":
	plt.title('Simulated Wright-Fisher Model: Recessive')
plt.legend(loc='lower right') #https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
plt.show()
#If the model parameter is "dominant", add one label, whereas if the model parameter is "recessive," add another label.
if model == "dominant":
		plt.savefig("allelic_frequency_vs_generations_dominant.png")
elif model == "recessive":
		plt.savefig("allelic_frequency_vs_generations_recessive.png")


