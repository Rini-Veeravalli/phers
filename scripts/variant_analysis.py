import pandas as pd
import numpy as np
from io import StringIO
import string


datalist = [[]]
#print(type(datalist))


#create empty dataframe
df_samplesGT = pd.DataFrame(columns=['chr3'])

with open("extracted_ClinVar_variants_Brugada.txt") as infile:
	count = 0
	for line in infile:
		#print(type(line))
		templist =[[]]
		
		if (count == 0):
			colnames = line.split("\t")
			count = count + 1
		else:
			rowentries = line.split("\t")
			#print(len(rowentries)-8)
			#print(len(datalist[0]))
			# TODO remove \n
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
			
			#print(rowentries[-1:])
			#print(templist[0:9])
			datalist.append(templist)
			


#print(colnames)
#print(datalist[0])
datalist.remove(datalist[0]) # removes empty [] at start of list

#print(len(colnames))
#print(len(datalist))
#print(len(datalist[0]))
length = len(datalist)
for i in range(length):
	datalist[i].remove(datalist[i][9]) # removes empty [] at end of row within item of list

#print(len(colnames))
#print(len(datalist))
#print(len(datalist[0]))
		
		#break
		
df_samplesGT = df_samplesGT.drop(columns=['chr3'])
df_finalsamplesGT = pd.DataFrame()
# remove \n in the last row of data
for col in df_samplesGT.iteritems():
	#print(col)
	#print(col[0])
	#print(col[1])
	#print(col[-1])
	#print((col[-1].map(lambda x : x.rstrip('\n'))))
	# add to new dataframe
	df_finalsamplesGT[str(col[0])] = col[-1].map(lambda x : x.rstrip('\n'))
	

#print(df_samplesGT)
print(df_finalsamplesGT.head())
#print(df_finalsamplesGT.size())
df_finalsamplesGT.to_csv("samplesGT_Brugada.csv", encoding = 'utf-8', index = False)

df = pd.DataFrame(datalist, columns = colnames)
print(df)

#df.to_csv("extracted_ClinVar_variants_Brugada.csv", encoding = 'utf-8', index = False) 
