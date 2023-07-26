import networkx as nx
import numpy as np
import pandas as pd
import sys
import tqdm

parameter = sys.argv[1]

# load the edge list
df = pd.read_csv("9606.protein.physical.links.v11.5.txt", sep=" ")

# set a cutoff for score
cutoff_score = 600
df = df[df.combined_score > cutoff_score]

ensembl_map = pd.read_csv("Homo_sapiens.GRCh38.110.entrez.tsv", sep="\t").set_index("protein_stable_id")["gene_stable_id"].to_dict()
del ensembl_map['-']

def try_ensembl_map(ensembl):
  try:
    return ensembl_map[ensembl]
  except KeyError:
    return np.nan

# map the ensemblle protein id to a an entrex id for each edge in the PPI
df["gene1"] = df.protein1.map(lambda x:try_ensembl_map(x.split(".")[1]))
df["gene2"] = df.protein2.map(lambda x:try_ensembl_map(x.split(".")[1]))

df = df.dropna()

# build a networkx graph
ppi = nx.Graph()
for i,v in df.iterrows():
  ppi.add_edge(v["gene1"],v["gene2"])

# list out the nodes in the network
nodelist = list(ppi.nodes())
# 
unnormalized_propagator = nx.to_scipy_sparse_array(ppi, nodelist = nodelist)

import scipy.sparse

# make a diagonal matrix to rescale the rows of the matrix
# making a random walk adjacency matrix
degree_array = [1/ppi.degree(i) for i in nodelist]
degree_normalization = scipy.sparse.diags(degree_array)

# apply it to the unnormalize propagator to make
# the one-step random walk propagator
propagator = degree_normalization @ unnormalized_propagator

# build up the random walk with restart operatorr
beta = int(parameter)/100
random_walker = (1-beta)*scipy.sparse.eye(len(nodelist))

cum_propagator = scipy.sparse.eye(len(nodelist))
for i in range(6):
    print(i)
    cum_propagator = propagator @ cum_propagator
    random_walker += (beta**i)*cum_propagator

mutation_unstacked = pd.read_csv("HMF_somatic_PASS_GENE_HIGH.txt", sep="\t", header=None)
mutation_unstacked[4] = 1
mutation = pd.pivot_table(mutation_unstacked, index=0, columns=3, values=4).T
mutation = mutation.reindex(nodelist).fillna(0.0)

smoothed_mutations = []
for i in tqdm.tqdm(mutation.columns):
    smoothed_mutations.append(random_walker * mutation.loc[:,i])
smoothed_mutation = pd.DataFrame(np.vstack(smoothed_mutations).T, index = mutation.index, columns=mutation.columns)


smoothed_mutation.to_csv("smoothed_mutation_"+parameter+".csv")
