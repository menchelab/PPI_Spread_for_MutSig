import networkx as nx
import numpy as np
import pandas as pd
import sys
import tqdm
import seaborn as sns

import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt

ppi = nx.read_gml("data/ppi.gml")
# list out the nodes in the network
nodelist = list(ppi.nodes())

mutation_unstacked = pd.read_csv("data/HMF_somatic_PASS_GENE_HIGH.txt", sep="\t", header=None)
mutation_unstacked[4] = 1
mutation = pd.pivot_table(mutation_unstacked, index=0, columns=3, values=4).T
mutation = mutation.reindex(nodelist).fillna(0.0)

# removed rows and columns with all zeros
# representing samples with no mutations in PPI genes and
# isolated, unmutated genes in the PPI respectively
mutation = mutation.T[(mutation**2).sum(axis=0) > 0].T
mutation = mutation[(mutation**2).sum(axis=1) > 0]


aaa = sns.clustermap(mutation.T, cmap="seismic", metric = 'cosine', col_cluster=False, center=0, vmax = 0.1, vmin=-0.1)
aaa.savefig("figures/original_clustermap.png")
aaa.fig.show()

import sklearn.decomposition
pca = sklearn.decomposition.PCA(n_components=100)
pca_transformed = pca.fit_transform(mutation.T)
pca_transformed = pd.DataFrame(pca_transformed, index=mutation.columns)



aaa = sns.clustermap(pca_transformed, cmap="seismic", metric = 'cosine', col_cluster=False, center=0, vmax = 0.1, vmin=-0.1)
aaa.savefig("figures/pca_original_clustermap.png")
aaa.fig.show()
