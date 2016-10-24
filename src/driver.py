
import re
import copy
import itertools
import networkx as nx

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
exit()

print("Calculate Degree Centrality...")
GeneDegrees = dict()

for gene_pair in PanCancerGenePairs:
    GeneDegrees[gene_pair[0]] = GeneDegrees.get(gene_pair[0],0) + 1
    GeneDegrees[gene_pair[1]] = GeneDegrees.get(gene_pair[1],0) + 1

y_star = max(GeneDegrees.values())
h_degree = sum([y_star - x for x in GeneDegrees.values()])
print("{0}".format(h_degree))

print("Calculate Closeness Centrality...")
def dist(source,dest):
    distance = 1;
    if (source,dest) in PanCancerGenePairs:
        return distance

    frontier = set([source])
    traversed = set()
    neighbors = set()

    while distance < DEPTH_LIM:
        distance += 1
        for cur_gene in frontier:
            for cur_pair in PanCancerGenePairs:
                if cur_pair[0] == cur_gene:
                    if cur_pair[1] == dest:
                        return distance
                    if not cur_pair[1] in traversed:
                        neighbors.add(cur_pair[1])
                        traversed.add(cur_pair[1])
                elif cur_pair[1] == cur_gene:
                    if cur_pair[0] == dest:
                        return distance
                    if not cur_pair[0] in traversed:
                        neighbors.add(cur_pair[0])
                        traversed.add(cur_pair[0])
        frontier = copy.deepcopy(neighbors)
        neighbors.clear()
    #print("depth exceeded limit ({0}) when finding path from {1} to {2}".format(distance,source,dest))
    return float("inf")

GeneDistances = dict()

for pair in itertools.product(GeneSubset, repeat=2):
    if pair[0] != pair[1]:
        GeneDistances[pair] = dist(pair[0],pair[1])

h_closeness = sum([1/x for x in GeneDistances.values()])
print("{0}".format(h_closeness))

print("Calculate Betweenness Centrality...")
def shortestPath(source,dest,path_len):
    frontier = set([source])
    traversed = set()
    neighbors = set()
    parent_maps = list()
    complete = False

    while not complete:
        parents = dict()
        for cur_gene in frontier:
            for cur_pair in PanCancerGenePairs:
                if cur_gene == cur_pair[0]:
                    parents[cur_pair[1]] = cur_gene
                    if cur_pair[1] == dest:
                        complete = True
                    elif cur_pair[1] not in traversed:
                        neighbors.add(cur_pair[1])
                        traversed.add(cur_pair[1])
                elif cur_gene == cur_pair[1]:
                    parents[cur_pair[0]] = cur_gene
                    if cur_pair[0] == dest:
                        complete = True
                    elif cur_pair[0] not in traversed:
                        neighbors.add(cur_pair[0])
                        traversed.add(cur_pair[0])

        parent_maps.append(copy.deepcopy(parents))
        frontier = copy.deepcopy(neighbors)
        neighbors.clear()

    sp = list()
    cur_sp_vertex = dest
    bottom_parent_idx = len(parent_maps) - 1
    while bottom_parent_idx >= 0:
        sp.append(cur_sp_vertex)
        #print("parent_maps[{0}]: {1}".format(bottom_parent_idx,parent_maps[bottom_parent_idx]))
        #print("parent_maps[{0}][{1}]: {2}".format(bottom_parent_idx,cur_sp_vertex,parent_maps[bottom_parent_idx][cur_sp_vertex]))
        cur_sp_vertex = parent_maps[bottom_parent_idx][cur_sp_vertex]
        bottom_parent_idx -= 1

    return sp

def btw(vertex):
    denominator = 0.0
    numerator = 0.0
    for pair in GeneDistances:
        path_len = GeneDistances[pair]
        if path_len < float("inf") and not vertex == pair[0] and not vertex == pair[1] and not pair[0] == pair[1]:
            denominator += 1
            if vertex in shortestPath(pair[0],pair[1],path_len):
                numerator += 1

    ratio = numerator/denominator
    #print("{0}:{1} = {2}".format(numerator,denominator,ratio))
    return(ratio)

GeneBetweenness = dict()

count = 0
for gene in GeneSubset:
    count += 1
    GeneBetweenness[gene] = btw(gene)
    print("GeneBetweeness[{0}] = {1} ({2} of {3})".format(gene,GeneBetweenness[gene],count,len(GeneSubset)))

h_betweenness = sum([x for x in GeneBetweenness.values()])
print("{0}".format(h_betweenness))

print("Calculate Eigenvector Centrality...")
def neighbors(gene):
    count = 0
    for pair in PanCancerGenePairs:
        if gene == pair[0] and not gene == pair[1]:
            count += 1
        elif gene == pair[1] and not gene == pair[0]:
            count +=1

    return count

GeneEigenvector = dict()

for gene in GeneSubset:
    GeneEigenvector[gene] = neighbors(gene)

h_eigenvector = sum([x for x in GeneEigenvector.values()])
print("{0}".format(h_eigenvector))