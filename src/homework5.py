# coding=utf-8
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

# return distance from root to leaf
def bfs(root,leaf):
    pass

# return number of components in graph
def num_components(graph):
    pass

# return diameter of sub-network
def diameter(sub_network):
    pass

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