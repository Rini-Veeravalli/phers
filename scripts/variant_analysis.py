# part I of variant analysis:
# script that transforms extracted_ClinVar_variants_[disease].txt:
#                        columns = Chrom, Pos, Ref, Alt, INFO_AC, INFO_AN, INFO_NS, INFO_AF, SAMPLE=GT (genetic carrier information of all tested participants (in one cell)
#                        rows = extracted ClinVar varianst of the disease
#  
#                   into 1) samplesGT_[disease].csv: columns = extracted ClinVar variants, rows = genetic variant results per participant 
#                        2) extracted_ClinVar_variants_[disease].csv: .csv version of original .txt file
# to be used to find carriers fur further analysis in variant.analysis.R

import pandas as pd
import numpy as np
from io import StringIO
import string


datalist = [[]]

#create empty dataframe
df_samplesGT = pd.DataFrame(columns=['chr3'])

with open("extracted_ClinVar_variants_Brugada.txt") as infile:
	count = 0
	for line in infile:
		templist =[[]]
		
		if (count == 0):
			colnames = line.split("\t")
			count = count + 1
		else:
			rowentries = line.split("\t")
			templist.insert(0, rowentries[0] )
			templist.insert(1, rowentries[1] )
			templist.insert(2, rowentries[2] )
			templist.insert(3, rowentries[3] )
			templist.insert(4, rowentries[4] )
			templist.insert(5, rowentries[5] )
			templist.insert(6, rowentries[6] )
			templist.insert(7, rowentries[7] )
			templist.insert(8, rowentries[-(len(rowentries)-8):])
			
			#add new col to samplesGT df with samples as rows
			df_samplesGT[str(rowentries[1])] = pd.Series(rowentries[-(len(rowentries)-8):])
			
			datalist.append(templist)

# dataframe to contain only genetically tested samples (row) for each variant (col)				
df_samplesGT = df_samplesGT.drop(columns=['chr3'])
df_finalsamplesGT = pd.DataFrame()

# remove \n in the last row of data
for col in df_samplesGT.iteritems():
	# add to new dataframe
	df_finalsamplesGT[str(col[0])] = col[-1].map(lambda x : x.rstrip('\n'))	
print(df_finalsamplesGT.head())

# save dataframe as .csv file, will be used to find carriers of each variant
df_finalsamplesGT.to_csv("samplesGT_Brugada.csv", encoding = 'utf-8', index = False)


# [optional] convert original extracted variant .txt file to .csv 
datalist.remove(datalist[0]) # removes empty [] at start of list
length = len(datalist)
for i in range(length):
	datalist[i].remove(datalist[i][9]) # removes empty [] at end of row within item of list
df = pd.DataFrame(datalist, columns = colnames)
print(df)

# save disease's extracted variants as .csv file
df.to_csv("extracted_ClinVar_variants_Brugada.csv", encoding = 'utf-8', index = False) 
