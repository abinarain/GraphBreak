#####################################################################################################################################################################
#Copyright PLease take permission before modifying or using this script
#Author Abhishek N. Singh
#Email abhishek.narain@iitdalumni.com
#Usage python3 programName inputFileName
#Date 19th June 2019
#Example python3 orderCommunity.py  communitiesCentral.csv
#Note the inputFileName should have a column by name node and another by name community separated by space or tab
#####################################################################################################################################################################
import pandas 
import numpy as np
import sys

df= pandas.read_csv(sys.argv[1], sep='\s+')
df = df.assign(G=df.groupby('community').community.transform('count')).sort_values(['G','community'],ascending=[False,True]).drop('G',1)
df.to_csv("communitiesOrdered.csv", sep=" ", index=False)
