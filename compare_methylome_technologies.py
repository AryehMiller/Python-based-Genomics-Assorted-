#!/usr/bin/env python3

"""
usage: python3 compare_methylome_technologies.py <MeDIP-seq RPKM bed> <MRE-seq RPKM bed> <WGBS methylation level bed>
Arguments:  <WGBS bed> WGBS_CGI_methylation.bed -- a BED file with the columns 1) chromosome, 2) start, 3) stop, 4) CGI name, and 5) average CpG methylation in CGI
output: several plots defining the correlation between two bedfiles and printed correlation coefficients
"""
import sys
import os
import matplotlib.pyplot as plt
from scipy import stats
###
##need 3 input files
if len(sys.argv) != 4:
	print(__doc__)
	sys.exit('error: incorrect number of input parameters, expecting four')

###
##Read in BED file
MeDIP = open(sys.argv[1],"r")
MeDIP_path = sys.argv[1]

# Creating dynamic output filenames
#Get the name
name_MeDIP = os.path.basename(MeDIP_path)
#Get the pruned name (w/o file ending)
shortname_MeDIP = os.path.splitext(name_MeDIP)[0]

##Read in BED file
MRE = open(sys.argv[2],"r")
MRE_path = sys.argv[2]

# Creating dynamic output filenames
#Get the name
name_MRE = os.path.basename(MRE_path)
#Get the pruned name (w/o file ending)
shortname_MRE = os.path.splitext(name_MRE)[0]

##Read in BED file
WGBS = open(sys.argv[3],"r")
WGBS_path = sys.argv[3]

# Creating dynamic output filenames
#Get the name
name_WGBS = os.path.basename(WGBS_path)
#Get the pruned name (w/o file ending)
shortname_WGBS = os.path.splitext(name_WGBS)[0]

#NOTE: These coordinates aren't mutual, so need to remove to allow for comparison
#15035943	15036236
#15038195	15038815
#15052974	15053360
#15068825	15070104
#15077051	15077824

# initialize list
Corr_Table_MeDIP = []

## SAME FOR LOOP STRUCTURE AS PREVIOUS SCRIPTS WHERE LINE IS STRIPPED AND APPENDED

#For loop to obtain 4th column
for line in MeDIP:
	each_line = line.strip().split()
	if each_line[1] not in ('15035943','15038195','15052974','15068825','15077051'):
		Corr_Table_MeDIP.append(float(each_line[4]))
		
# initialize list
Corr_Table_MRE = []

## SAME FOR LOOP STRUCTURE AS PREVIOUS SCRIPTS WHERE LINE IS STRIPPED AND APPENDED

#For loop to obtain 4th column
for line in MRE:
	each_line = line.strip().split()
	if each_line[1] not in ('15035943','15038195','15052974','15068825','15077051'):
		Corr_Table_MRE.append(float(each_line[4]))


## SAME FOR LOOP STRUCTURE AS PREVIOUS SCRIPTS WHERE LINE IS STRIPPED AND APPENDED

#For loop to obtain 4th column
# initialize list
Corr_Table_WGBS = []
for line in WGBS:
	each_line = line.strip().split()
	Corr_Table_WGBS.append(float(each_line[5]))
	
#plot source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html

##SAME ANNOTATION FOR ALL PLOTS

## clear new plot
## plot x and y in lightskyblue at alpha 0.5
#Label x axis
#Label y axis
#Label main title
#Save fig
#Print the spearman correlation coefficient


## MRE vs MeDIP
fig = plt.figure() 
plt.scatter(Corr_Table_MRE, Corr_Table_MeDIP, color = "lightskyblue", alpha=0.5)
plt.title('MRE vs. MeDIP Methylation Scores')
plt.xlabel('MRE Methylation Scores')
plt.ylabel('MeDIP Methylation Scores')
plt.savefig(shortname_MRE + "_vs_" + shortname_MeDIP + ".png")
print("Correlation between MRE and MeDIP:")
print(stats.spearmanr(Corr_Table_MRE, Corr_Table_MeDIP))

## MRE vs WGBS
fig = plt.figure()
plt.scatter(Corr_Table_WGBS, Corr_Table_MRE, color = "lightskyblue", alpha = 0.5)
plt.title('MRE vs. WGBS Methylation Scores')
plt.xlabel('WGBS Methylation Scores')
plt.ylabel('MRE Methylation Scores')
fig.savefig(shortname_MRE + "_vs_" + shortname_WGBS + ".png")
print("Correlation between WGBS and MRE:")
print(stats.spearmanr(Corr_Table_WGBS, Corr_Table_MRE))

## MeDIP vs WGBS
fig = plt.figure()
plt.scatter(Corr_Table_WGBS, Corr_Table_MeDIP, color = "lightskyblue", alpha=0.5)
plt.title('MeDIP vs. WGBS Methylation Scores')
plt.xlabel('WGBS Methylation Scores')
plt.ylabel('MeDIP Methylation Scores')
plt.savefig(shortname_MeDIP + "_vs_" + shortname_WGBS + ".png")
print("Correlation between WGBS and MeDIP:")
print(stats.spearmanr(Corr_Table_WGBS, Corr_Table_MeDIP))

##Outliers
# Find outliers-- these will be them
#print(max(Corr_Table_MeDIP))
#print(max(Corr_Table_MRE))
#WGBS 1.0 point probably not an outlier
#984.6279
#66384.4794
#Searched the file to find the coordinates of those outliers-- they are below-- now need to remove for no-outlier comparison
#chr21	9825442	9826296	305.6206	984.6279
#chr21	9825442	9826296	11779.8595	66384.4794

#NO OUTLIERS

### NOTE: I re-paste this part because to overwrite the previous portion and allow for the current coding scheme below. Should be OK, although I'm sure there's a far easier way.

##Read in BED file
MeDIP = open(sys.argv[1],"r")
MeDIP_path = sys.argv[1]

# Creating dynamic output filenames
#Get the name
name_MeDIP = os.path.basename(MeDIP_path)
#Get the pruned name (w/o file ending)
shortname_MeDIP = os.path.splitext(name_MeDIP)[0]

##Read in BED file
MRE = open(sys.argv[2],"r")
MRE_path = sys.argv[2]

# Creating dynamic output filenames
#Get the name
name_MRE = os.path.basename(MRE_path)
#Get the pruned name (w/o file ending)
shortname_MRE = os.path.splitext(name_MRE)[0]

##Read in BED file
WGBS = open(sys.argv[3],"r")
WGBS_path = sys.argv[3]

# Creating dynamic output filenames
#Get the name
name_WGBS = os.path.basename(WGBS_path)
#Get the pruned name (w/o file ending)
shortname_WGBS = os.path.splitext(name_WGBS)[0]


# initialize list
Corr_Table_MeDIP_NoOutliers = []

#For loop to obtain 4th column
for line in MeDIP:
	each_line = line.strip().split()
	if each_line[1] not in ('9825442','15035943','15038195','15052974','15068825','15077051'):
		Corr_Table_MeDIP_NoOutliers.append(float(each_line[4]))

# initialize list
Corr_Table_MRE_NoOutliers = []

#For loop to obtain 4th column
for line in MRE:
	each_line = line.strip().split()
	if each_line[1] not in ('9825442','15035943','15038195','15052974','15068825','15077051'):
		Corr_Table_MRE_NoOutliers.append(float(each_line[4]))

# initialize list
Corr_Table_WGBS_NoOutliers = []

#For loop to obtain 4th column
for line in WGBS:
	each_line = line.strip().split()
	if each_line[1] not in ('9825442'):
		Corr_Table_WGBS_NoOutliers.append(float(each_line[5]))

##SAME ANNOTATION FOR ALL PLOTS

## clear new plot
## plot x and y in lightskyblue at alpha 0.5
#Label x axis
#Label y axis
#Label main title
#Save fig
#Print the spearman correlation coefficient

#plot source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
## MRE vs MeDIP	
fig = plt.figure() #start new fig
plt.scatter(Corr_Table_MRE_NoOutliers, Corr_Table_MeDIP_NoOutliers, color = "lightskyblue", alpha=0.5) #x y
plt.title('MRE vs. MeDIP Methylation Scores Outliers Removed')
plt.xlabel('MRE Methylation Scores')
plt.ylabel('MeDIP Methylation Scores')
plt.savefig(shortname_MRE + "_vs_" + shortname_MeDIP + "_outliers_removed.png")
print("Correlation between MRE and MeDIP Outliers Removed:")
print(stats.spearmanr(Corr_Table_MRE_NoOutliers, Corr_Table_MeDIP_NoOutliers))

## MRE vs WGBS
fig = plt.figure()
plt.scatter(Corr_Table_WGBS_NoOutliers, Corr_Table_MRE_NoOutliers, color = "lightskyblue", alpha = 0.5)
plt.title('MRE vs. WGBS Methylation Scores Outliers Removed')
plt.xlabel('WGBS Methylation Scores')
plt.ylabel('MRE Methylation Scores')
fig.savefig(shortname_MRE + "_vs_" + shortname_WGBS + "_outliers_removed.png")
print("Correlation between WGBS and MRE Outliers Removed:")
print(stats.spearmanr(Corr_Table_WGBS_NoOutliers, Corr_Table_MRE_NoOutliers))

## MeDIP vs WGBS
fig = plt.figure()
plt.scatter(Corr_Table_WGBS_NoOutliers, Corr_Table_MeDIP_NoOutliers, color = "lightskyblue", alpha=0.5)
plt.title('MeDIP vs. WGBS Methylation Scores Outliers Removed')
plt.xlabel('WGBS Methylation Scores')
plt.ylabel('MeDIP Methylation Scores')
plt.savefig(shortname_MeDIP + "_vs_" + shortname_WGBS + "_outliers_removed.png")
print("Correlation between WGBS and MeDIP Outliers Removed:")
print(stats.spearmanr(Corr_Table_WGBS_NoOutliers, Corr_Table_MeDIP_NoOutliers))

