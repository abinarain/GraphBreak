# ReGen
## Regulatory Gene Network Analysis Toolkit 

Download via Git: git clone https://github.com/abinarain/ReGen.git

#CopyRight: This software tool is a copyright of the author. Please take permission before use or modification. License of use is as per MIT License agreement.

#Author: Abhishek Narain Singh

#Email: abhishek.narain@iitdalumni.com

#Organization affiliation: A.I. Virtanen Institute for Molecular Sciences, University of Eastern Finland, Kuopio, Finland

Dependencies: OS - Linux, Python 3.7.3, OpenMP, Additional Python Modules needed to be installed: networkx, matplotlib and networkit

Install cython, sklearn, seaborn and cmake which are Networkit dependencies:
sudo apt install cmake
pip3 install cython 
pip3 install seaborn
pip3 install sklearn

Then to install networkit:
pip3 install --user networkit

Alternatively, this can be built from source, for which you can 1st download it from GitHub:
git clone https://github.com/networkit/networkit.git
cd networkit
sudo python3 setup.py build_ext
#Note here again cmake and cython would be needed for installation


To install networkx: pip3 install --user networkx

To install matplotlib: pip3 install --user matplotlib

To load Python 3.7.3 using module: module load python/3.7.3


## For Louvaine Algorithm Parallel implementation detecting communities:

Sample data provided for the purpose small1976.txt (this has been picked up from R database), GTExHead2000Artery_Aorta.v7.signif_variant_gene_pairs.txt (This is the top 2000 genomic variants to Gene expression significant association from the GTEx V7 file that is publicly available to download from their portal).

#Usage: python programName NumberbelowwhichToColorDifferently spaceDelimitedDataFileWithColumnNames column1Name column2Name NumberOfParallelCoresOpenMP

 Example: python3 bipartiteLinearKit.py 2 small1976.txt pollinator plant 2
          #For generating cluster without any bipartite plot (to save tremendous amount of time)
          python3 bipartiteLinearKitNoAbi.py 2 small1976.txt pollinator plant 2
 
 Example: python bipartiteLinearKit.py 2 GTExHead2000Artery_Aorta.v7.signif_variant_gene_pairs.txt variant_id gene_id 5

The output file generated would lead to this kind of verbose print:

Created the Graph Structure

Colored the Nodes

Drawing prepared

Graph Plotted by name abiPlot.png

PLM(balanced,pc,turbo) detected communities in 0.002198457717895508 [s]

solution properties:

-------------------  ----------

#communities         15

min community size     2

max community size   283

avg. community size   81.8667

modularity             0.804935

-------------------  ----------

0.8049354842503793

wrote communities to: ./communities.partition

Communities Plotted by name communityPlot.png

Prepared the Communities file as nodesCommunity.txt

The output files will comprise of:

 nodesCommunity.txt* - This comprises of a list with all the entries assigned to different communities and their numbering
 
 communityPlot.png* - plot of various communities
 
 communities.partition* 
 
 abiPlot.png* - Plot in Bipartite Linear fashion with color difference the first argument after the program name
 
 nodes.txt* - list of nodes

## For Girvan-Newman Algorithm for Community Detection:

Note: This is not a parallel implementation, and so should not be preferred for large datasets. 

#Usage: python programName degreeAboveOrEqualToWhichToColorDifferentlyIn fileName 1stColumnName 2ndColumnName

 Example: python3 bipartiteGN.py 2 small1976.txt pollinator plant 
 
 Example: python3 bipartiteGN.py 2 GTExHead2000Artery_Aorta.v7.signif_variant_gene_pairs.txt variant_id gene_id

The output would look like:

Created the Graph

Colored the Nodes

Drawing for Plot prepared

Graph Plotted by name abiPlot.png

The list of communities are successfully fed in communitiesCentral.csv file

Output files generated:

communitiesCentral.csv - list of nodes with the corresponding group number for the community each one belongs to

abiPlot.png - plot of the bipartite nature of the connectedness of graph nodes


## Extract Community Members (Genes/Transcript or Variants - SNPs) with a Prefix

Usage: python3 extractPrefix.py nodesCommunity.txt Prefix

#Note the inputFileName nodesCommunity.txt should have a column by name node and another by name community separated by space or tab
#The output file is by name communitiesPrefix.csv

Example to get all the list of genes in various communities:
python3 extractPrefix.py nodeCommunity.txt ENS

## Order Community Members (Genes/Transcript or Variants - SNPs) in Descending order

python3 orderCommunity.py communitiesPrefix.csv #The columns should be named as node and community separated by tab or space


## Then the downstream analysis is by means of www.reactome.org

website of Author: www.tinyurl.com/abinarain 
