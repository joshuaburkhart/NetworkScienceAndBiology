
import re
import networkx as nx
import json

EXTRACT_UIDS_RGX = '^[^\t]+\t[^\t]+\tuniprotkb:[^_]+_(?P<UIDA>[^\(]+)\(shortlabel\)\tuniprotkb:[^_]+_(?P<UIDB>[^\(]+)\(shortlabel\)\t'
EXTRACT_GENE_RGX = '^(?P<UID>.+)$'
INTERACTIONS_FILE_PATH = '/Users/joshuaburkhart/SoftwareProjects/NetworkScienceAndBiology/input/homo_sapiens.mitab.interactions.txt'
GENE_SUBSET_FILE_PATH = '/Users/joshuaburkhart/SoftwareProjects/NetworkScienceAndBiology/input/TCGA_PanCancer_Nature_ 24132290.txt'
DEPTH_LIM = 25

print("Extract interactions into pair-wise gene symbols...")
GeneAdjMatrix = set()

in_fptr = open(INTERACTIONS_FILE_PATH)
while 1:
    line = in_fptr.readline()
    if not line:
        break
    match = re.match(EXTRACT_UIDS_RGX, line)
    if match:
        GeneAdjMatrix.add((match.group('UIDA'), match.group('UIDB')))
in_fptr.close()

print("Map the TCGA PanCancer genes...")
GeneSubset = set()

in_fptr = open(GENE_SUBSET_FILE_PATH)
while 1:
    line = in_fptr.readline()
    if not line:
        break
    match = re.match(EXTRACT_GENE_RGX,line)
    if match:
        GeneSubset.add(match.group('UID'))
in_fptr.close()

PanCancerGenePairs = [gene_pair for gene_pair in GeneAdjMatrix if gene_pair[0] in GeneSubset and gene_pair[1] in GeneSubset]

G = nx.Graph()
for pair in PanCancerGenePairs:
    G.add_edge(pair[0],pair[1])

nx.write_graphml(G,"/Users/joshuaburkhart/tmp/G.xml")

print("Calculate Degree Centrality...")
output = open("/Users/joshuaburkhart/tmp/degree_centrality_G.json", 'w')
json.dump(nx.degree_centrality(G), output)
output.close()

print("Calculate Eigenvector Centrality...")
output = open("/Users/joshuaburkhart/tmp/eigenvector_centrality_G.json", 'w')
json.dump(nx.eigenvector_centrality(G), output)
output.close()

print("Calculate Katz Centrality...")
output = open("/Users/joshuaburkhart/tmp/katz_centrality_G.json", 'w')
json.dump(nx.katz_centrality(G,alpha=0.03,max_iter=10000), output)
output.close()

print("Calculate PageRank Centrality...")
output = open("/Users/joshuaburkhart/tmp/pagerank_centrality_G.json", 'w')
json.dump(nx.pagerank(G), output)
output.close()

print("Done.")
