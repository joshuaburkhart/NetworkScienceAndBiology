# coding=utf-8
import re
import math
import random
import networkx as nx
from collections import deque

EXTRACT_UIDS_RGX = '^[^\t]+\t[^\t]+\tuniprotkb:[^_]+_(?P<UIDA>[^\(]+)\(shortlabel\)\tuniprotkb:[^_]+_(?P<UIDB>[^\(]+)\(shortlabel\)\t'
EXTRACT_GENE_RGX = '^(?P<UID>.+)$'
INTERACTIONS_FILE_PATH = '/home/burkhart/Software/NetworkScienceAndBiology/input/homo_sapiens.mitab.interactions.txt'
GENE_SUBSET_FILE_PATH = '/home/burkhart/Software/NetworkScienceAndBiology/input/TCGA_PanCancer_Nature_ 24132290.txt'

random.seed(88)

print("Extract interactions into pair-wise gene symbols...")
GeneAdjMatrix = set()

# create circle network
def circle_network(n,c):
    circle = nx.Graph()
    for i in range(n):
        if i % (c + 1) >= (c + 1) // 2:
            k = (i + ((c + 1) // 2)) % n
            print('+ special edge: {0} -> {1}'.format(i,k))
            circle.add_edge(i,k)
        for j in range(c - 1):
            k = (i + (j // 2 + 1) * pow(-1, j % 2)) % n
            print('+ edge: {0} -> {1}'.format(i,k))
            circle.add_edge(i,k)
    return circle

# create permuted circle network
def permuted_circle_network(n,c,p):
    num_permutations = math.floor(n * p)
    permuted_circle = circle_network(n,c)
    for i in range(num_permutations):
        rand_edge = random.choice(permuted_circle.edges())
        print('- edge: {0} -> {1}'.format(rand_edge[0],rand_edge[1]))
        permuted_circle.remove_edge(rand_edge[0],rand_edge[1])
        j_add = 0
        k_add = 0
        while k_add == j_add or (j_add,k_add) in permuted_circle.edges():
            j_add = random.randint(0,n)
            k_add = random.randint(0,n)
        print('-/+ edge: {0} -> {1}'.format(j_add,k_add))
        permuted_circle.add_edge(j_add,k_add)

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

circle0 = circle_network(20,6)
nx.write_graphml(circle0,"/home/burkhart/Software/NetworkScienceAndBiology/output/circle.xml")

pcircle0 = permuted_circle_network(20,6,0.2)
nx.write_graphml(circle0,"/home/burkhart/Software/NetworkScienceAndBiology/output/pcircle.xml")

