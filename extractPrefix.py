#####################################################################################################################################################################
#Copyright PLease take permission before modifying or using this script
#Author Abhishek N. Singh
#Description: This script extracts all the rows match a certain prefix in node column
#Email abhishek.narain@iitdalumni.com
#Usage python3 programName inputFileName prefix
#Date 27th June 2019
#-bash-4.2$ python3 extractPrefix.py communitiesOrderedSignificantArteryAortaLiteGTExv7.csv ENS
#Example python3 extractPrefix.py  communitiesCentral.csv ENS
#Note the inputFileName should have a column by name node and another by name community separated by space or tab
#The output file is by name communitiesPrefix.csv
#####################################################################################################################################################################
import pandas
import numpy as np
import sys
import re

df= pandas.read_csv(sys.argv[1], sep='\s+') #Here we read the filename into dataframe
prefix = sys.argv[2] #Storing the prefix variable
#print(df.node)
#print(re.findall(r"^prefix",df.node))
mask = df['node'].str.match(prefix)
#print(df[mask])
df[mask].to_csv("communitiesPrefix.csv", sep=" ", index=False)
#df.to_csv("communitiesOrdered.csv", sep=" ", index=False)
