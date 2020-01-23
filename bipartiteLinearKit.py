#######################################################################################################################################################
#CopyRight: This software tool is a copyright of the author. Please take permission before use or modification.
#Author: Abhishek Narain Singh
#Description: This code is for Bipartite graph plotting. For example SNPs and different Phenotypes, or SNPs and different Genes, such as in eQTL
#Email: abhishek.narain@iitdalumni.com
#Example: -bash-4.2$ python3 bipartiteLinearKit.py 2 small1976.txt pollinator plant 2
# Usage: python3 programName NumberbelowwhichToColorDifferently spaceDelimitedDataFileWithColumnNames column1Name column2Name NumberOfParallelCoresOpenMP
#Example: python3 bipartiteLinearKit.py 2 ~/GTEx/GTEx_Analysis_v7_eQTL/Artery_Aorta.v7.signif_variant_gene_pairs.txt variant_id gene_id 30
#Date: 25th June 2019
########################################################################################################################################################
import matplotlib
#matplotlib.use('QT4Agg')
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import sys
import pandas as pd
from networkx.algorithms import community
import networkit as nk
import subprocess

df = pd.read_csv(sys.argv[2],sep='\s+') #Here goes the file name which is space or tab separated and 1st row as names of the columns
item1 = df[sys.argv[3]].unique() #Here goes the first column variable name such as the name of the genes
#print(item1)
item2 = df[sys.argv[4]].unique() #Here goes the second column variable name such as the name of the SNPs
#print(item2)
edges1 = df[sys.argv[3]]
edges2 = df[sys.argv[4]]
edges = pd.concat([edges1, edges2], axis=1)
edgesArray = edges.values

B = nx.Graph()
#SNPs = [1,2,3,4]
#Genes = ['a','b','c']
#Edge_Weight = ('r','r','b','b','g','r') #These weight are -log base 10 of the p-value for association of a SNP to Gene
B.add_nodes_from(item1, bipartite=0) # Add the node attribute "bipartite"
B.add_nodes_from(item2, bipartite=1)
#edges = [[1,'a'], [1,'b'], [2,'b'], [2,'c'], [3,'c'], [4,'a']]
B.add_edges_from(edgesArray)

print("Created the Graph Structure")
#Separating the nodes by group
r = {n for n, d in B.nodes(data=True) if d['bipartite']==0}  #Getting the top nodes
l = set(B) - r #Getting the lower nodes

#print(set(B))#Sets store unordered values so this is not needed
#Creating a File where the list of Nodes are written
listOfNodes = list(B.nodes)
with open ("nodes.txt", "w") as thefile:
    for item in listOfNodes:
        thefile.write(item + "\n")
thefile.close()
#print(list(B.nodes)[3])
#print(list(B))
#print(list(B.nodes(data=True)))
# Separate by group
#l, r = nx.bipartite.sets(B)
pos = {}

#print(l)
#print(r)
# Update position for node from each group THis will be needed for two parallel lines as bipartite
pos.update((node, (1, index)) for index, node in enumerate(l))
#print(pos)
pos.update((node, (2, index)) for index, node in enumerate(r))
#print(pos)
color_map = []
for nodeCount in range(len(item1)):
    #print(nodeCount)
    color_map.append('pink') #Item 1 objects colored one color
for nodeCount in range(len(item2)):
    #print(nodeCount)
    color_map.append('green') #Item 2 objects colored second color
print("Colored the Nodes")
#print(color_map)
#This is for two parallel line bipartite graph. Put pos=pos as an argument and see . To plot based on some edge weights
#nx.draw(B, pos=pos, with_labels=True, edge_color=Edge_Weight, node_color=color_map, node_size=1500, font_size=25, font_color="yellow", font_weight="bold",edge_cmap=plt.get_cmap('BuGn'), label ="SNP To Gene eQTL Associations Cis & Trans")
#To plot with degrees of association in linear bipartite
nx.draw(B, pos=pos, with_labels=True, edge_color=['blue' if B.degree[e[0]] >= int(sys.argv[1]) else 'red' for e in B.edges],font_size=4,font_weight="bold", node_color=color_map, font_color="black", edge_cmap=plt.get_cmap('BuGn'), label ="SNP To Gene eQTL Associations Cis & Trans")
#We make circular network plot with argument in command line for the degree of connectedness and above that needs to be colored differently
#nx.draw_circular(B,with_labels=True, edge_color=['blue' if B.degree[e[0]] >= int(sys.argv[1]) else 'red' for e in B.edges], node_color=color_map, font_color="black",edge_cmap=plt.get_cmap('Blues'), label ="SNP To Gene eQTL Associations Cis & Trans" )
#plt.title("SNP to Gene eQTL Association")
plt.title('ReGen Bipartite Plot', color='magenta')
#plt.show()
print("Drawing for Circular Plot prepared")
plt.savefig('abiPlot.png', bbox_inches='tight')
print("Graph Plotted by name abiPlot.png")
#plt.show()
plt.clf() #Clear the figure
####################################Community Detection is Done by using Network Kit Parallel Louvain's algorithm####################################
nk.setNumberOfThreads(int(sys.argv[5])) # Setting the number of Parallel Threads in OpenMP
nkG = nk.nxadapter.nx2nk(B, weightAttr=None) #Now nkG is the converted graph in networkit format
communities = nk.community.detectCommunities(nkG)

#nxG = nk.nxadapter.nk2nx(communities)

print(nk.community.Modularity().getQuality(communities, nkG))

#Write the community partitioning
nk.community.writeCommunities(communities, "./communities.partition")
#Plotting Communities (uncomment this if you want it)
nk.viztasks.drawCommunityGraph(nkG,communities)
plt.savefig('communityPlot.png', bbox_inches='tight')
print("Communities Plotted by name communityPlot.png")
#plt.show()
with open('nodesCommunity.txt', "w") as outfile:
    subprocess.call(["paste", "nodes.txt", "communities.partition"], stdout=outfile)
outfile.close()
#subprocess.call(["paste", "nodes.txt", "communities.partition", ">", "nodesCommunity.txt"], shell=True)
print("Prepared the Communities file as nodesCommunity.txt")
#c2 = communities.getMembers(3)
#print(c2)
#print(nk.getCurrentNumberOfThreads())
