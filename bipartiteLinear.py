#######################################################################################################################################################
#CopyRight: This software tool is a copyright of the author. Please take permission before use or modification.
#Author: Abhishek Narain Singh
#Description: This code is for Bipartite graph plotting. For example SNPs and different Phenotypes, or SNPs and different Genes, such as in eQTL
#Email: abhishek.narain@iitdalumni.com
#Usage: python programName degreeAboveOrEqualToWhichToColorDifferentlyIn Blue, lower than this gets Red
#Example: python3 bipartiteLinear.py 3 ~/GTEx/GTEx_Analysis_v7_eQTL/Artery_Aorta.v7.signif_variant_gene_pairs.txt variant_id gene_id
#Date: 12th June 2019
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

df = pd.read_csv(sys.argv[2],sep='\s+') #Here goes the file name which is space or tab separated and 1st row as names of the columns
item1 = df[sys.argv[3]].unique() #Here goes the first column variable name such as the name of the genes
item2 = df[sys.argv[4]].unique() #Here goes the second column variable name such as the name of the SNPs
edges1 = df[sys.argv[3]]
edges2 = df[sys.argv[4]]
edges = pd.concat([edges1, edges2], axis=1)
edgesArray = edges.values

B = nk.Graph()
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

#######################################################Here we Generate the Communities########################################
communities_generator = community.girvan_newman(B) #Finds communities in a graph using the Girvan–Newman Division method for centrality in packing
#communities_generator = community.greedy_modularity_communities(B)#Find communities in graph using Clauset-Newman-Moore greedy modularity maximization.
#communities_generator = community.asyn_fluidc(B) #Returns communities in G as detected by Fluid Communities algorithm.
#communities_generator = community.label_propagation_communities(B) #Generates community sets determined by label propagation
#communities_generator = community.kernighan_lin_bisection(B) #Partition a graph into two blocks using the Kernighan–Lin algorithm.
top_level_communities = next(communities_generator)
next_level_communities = next(communities_generator)

data = [[element, "GROUP-{}".format(ii + 1)] for ii, st in enumerate(next_level_communities) for element in sorted(st)]
#print(data)

frame = pd.DataFrame(data=data, columns=['node', 'community'])
frame.to_csv("communitiesCentral.csv", sep=" ", index=False)
print("The list of communities are successfully fed in communitiesCentral.csv file")
###################################Here we plot the community data#################################
#df = pd.read_csv("communitiesCentral.csv", sep='\s+')

#item1 = df['node'].unique()
#item2 = df.community.unique()
#edges1 = df.node
#edges2 = df.community
#edges = pd.concat([edges1, edges2], axis=1)
#edgesArray = edges.values


#B = nx.Graph()
#B.add_nodes_from(item1, bipartite=0) # Add the node attribute "bipartite"
#B.add_nodes_from(item2, bipartite=1)

#B.add_edges_from(edgesArray)

#nx.draw(B, with_labels=True, edge_color='green', node_color='pink', font_color="black", edge_cmap=plt.get_cmap('BuGn'), label ="SNP and Gene Communities")

#plt.title('Communities Plot')
#plt.show()
####################################################################################################

