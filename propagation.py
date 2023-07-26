import networkx as nx
import numpy as np
import pandas as pd
import sys
import tqdm

parameter = sys.argv[1]

ppi = nx.read_gml("data/ppi.gml")

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

random_walker = scipy.sparse.eye(len(nodelist))
cum_propagator = scipy.sparse.eye(len(nodelist))
for i in range(6):
    print(i)
    cum_propagator = propagator @ cum_propagator
    random_walker += (beta**i)*cum_propagator
random_walker = (1-beta)*random_walker

mutation_unstacked = pd.read_csv("data/HMF_somatic_PASS_GENE_HIGH.txt", sep="\t", header=None)
mutation_unstacked[4] = 1
mutation = pd.pivot_table(mutation_unstacked, index=0, columns=3, values=4).T
mutation = mutation.reindex(nodelist).fillna(0.0)

smoothed_mutations = []
for i in tqdm.tqdm(mutation.columns):
    smoothed_mutations.append(random_walker * mutation.loc[:,i])
smoothed_mutation = pd.DataFrame(np.vstack(smoothed_mutations).T, index = mutation.index, columns=mutation.columns)

# removed rows and columns with all zeros
# representing samples with no mutations in PPI genes and
# isolated, unmutated genes in the PPI respectively
smoothed_mutation = smoothed_mutation.T[(smoothed_mutation**2).sum(axis=0) > 0].T
smoothed_mutation = smoothed_mutation[(smoothed_mutation**2).sum(axis=1) > 0]

import sklearn.decomposition
# PCA to find the characteristic dimensions of
# "mutation space" 
pca = sklearn.decomposition.PCA(n_components=50)

# transfomred features are per sample and are the
# condensed representation of the mutational profile of each sample
transformed_features = pd.DataFrame(pca.fit_transform(smoothed_mutation), index=smoothed_mutation.index)
transformed_features.to_csv(f"data/pca_smooth_features_"+parameter+".csv")

# the compnents represent which genes'mutationns comttribute the 
# prrofile of eacch sample
components = pd.DataFrame(pca.components_, columns=smoothed_mutation.columns)
components.to_csv("data/pca_components_"+parameter+".csv")

draw=False
if draw:
    os.makedirs("figures/smoothing_"+parameter, exist_ok = True)

    # replace the ensembl ids in the collumnns with hgnc gene symbols
    smoothed_mutation_symbol = smoothed_mutation.copy()
    smoothed_mutation_symbol.columns = smoothed_mutation.columns.map(lambda x : ensembl_symbol.get(x, "Pass"))

    # plot the relatiosnhip between genes and samples after smoothng
    sns.clustermap(smoothed_mutation_symbol, vmax=0.05, metric='cosine')
    plt.savefig("figures/smoothing_"+parameter+"/clustermap_smoothed_mutation.png", dpi=500)
    plt.savefig("figures/smoothing_"+parameter+"/clustermap_smoothed_mutation.svg")
