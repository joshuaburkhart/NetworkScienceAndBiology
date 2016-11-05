# coding=utf-8
import re
import json
import networkx as nx
from collections import deque

EXTRACT_UIDS_RGX = '^[^\t]+\t[^\t]+\tuniprotkb:[^_]+_(?P<UIDA>[^\(]+)\(shortlabel\)\tuniprotkb:[^_]+_(?P<UIDB>[^\(]+)\(shortlabel\)\t'
EXTRACT_GENE_RGX = '^(?P<UID>.+)$'
INTERACTIONS_FILE_PATH = '/home/burkhart/Software/NetworkScienceAndBiology/input/homo_sapiens.mitab.interactions.txt'
GENE_SUBSET_FILE_PATH = '/home/burkhart/Software/NetworkScienceAndBiology/input/TCGA_PanCancer_Nature_ 24132290.txt'

print("Extract interactions into pair-wise gene symbols...")
GeneAdjMatrix = set()

# return dict of distances from root to leaf
def bfs(graph,vertex):
    depth = 0
    vertex_distances = dict()
    vertex_distances[vertex] = depth
    frontier = deque([vertex])

    while len(frontier) > 0:
        parent = frontier.popleft()
        depth = vertex_distances.get(parent)
        children = graph.neighbors(parent)
        for child in children:
            if vertex_distances.get(child) is None:
                vertex_distances[child] = (depth + 1)
                frontier.append(child)
    return vertex_distances

# return components in graph
def components(distances):
    unique_components = list()
    seen = set()
    for skey in distances.keys():
        if skey not in seen:
            new_component = list(distances[skey].keys())
            unique_components.append(new_component)
            seen.update(new_component)
    return unique_components

# return diameter of sub-network
def diameter(sub_network,distances):
    diam = 0
    for vertex in sub_network:
        max_dist = max(distances.get(vertex).values())
        diam = max_dist if max_dist > diam else diam
    return diam

def avg_s_pth(sub_network,distances):
    s_pth_dists = list()
    for vertex in sub_network:
        s_pth_dists.extend(list(distances.get(vertex).values()))
    return sum(s_pth_dists)/len(s_pth_dists)

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

dists = dict()

for vertex in G.nodes():
    dists[vertex] = bfs(G,vertex)

print("Writing distances...")
output = open('/home/burkhart/Software/NetworkScienceAndBiology/output/distances.txt', 'w')
for skey in sorted(list(dists.keys())):
    output.write('Distances from {0}\n'.format(skey))
    for tkey in sorted(list(dists[skey].keys())):
        output.write('\t{0}: {1}\n'.format(tkey,dists[skey].get(tkey,"Not Reachable")))
    output.write('\n')
output.close()

comps = components(dists)

print("Writing components...")
output = open('/home/burkhart/Software/NetworkScienceAndBiology/output/components.txt', 'w')
output.write('Discovered {0} components\n\n'.format(len(comps)))
for comp in sorted(list(comps),key=len):
    output.write('Unique component of size {0}\n'.format(len(comp)))
    output.write(' diameter {0}\n'.format(diameter(comp,dists)))
    output.write(' avg shortest path {0}\n'.format(avg_s_pth(comp,dists)))
    for element in sorted(list(comp)):
        output.write('\t{0}\n'.format(element))
    output.write('\n')
output.close()
