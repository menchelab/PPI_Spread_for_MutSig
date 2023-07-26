import pandas as pd
import json
import sklearn.decomposition as sd
import numpy as np

import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt

import os
import sys
parameter =  sys.argv[1]

# helper dictionary to map ensembl to hgnc symbol
name_mapper = json.load(open("data/protein-coding_gene.json"))['response']['docs']
ensembl_symbol = {i['ensembl_gene_id']:i['symbol'] for i in name_mapper if 'ensembl_gene_id' in i.keys()}

# removed rows and columns with all zeros
# representing samples with no mutations in PPI genes and
# isolated, unmutated genes in the PPI respectively
smoothed_mutation = pd.read_csv("data/smoothed_mutation_"+parameter+".csv", index_col=0).T 
smoothed_mutation = smoothed_mutation.T[(smoothed_mutation**2).sum(axis=0) > 0].T
smoothed_mutation = smoothed_mutation[(smoothed_mutation**2).sum(axis=1) > 0]

# PCA to find the characteristic dimensions of
# "mutation space" 
pca = sd.PCA(n_components=50)
transformed_features = pca.fit_transform(smoothed_mutation)

# 
os.makedirs("figures/smoothing_"+parameter, exist_ok = True) 

os.makedirs("figures/smoothing_"+parameter+"/pca/", exist_ok = True) 
i = 0
j = 1
fig,ax = plt.subplots()
ax.scatter(transformed_features[:,i], transformed_features[:,j])
ax.set_xlabel("PCA-"+str(j))
ax.set_ylabel("PCA-"+str(i))
fig.savefig("figures/smoothing_"+parameter+"/pca/PCA_"+str(i)+"_"+str(j)+".svg")
fig.savefig("figures/smoothing_"+parameter+"/pca/PCA_"+str(i)+"_"+str(j)+".png", dpi=500)

# load up the mutational signatures
signatures = pd.read_csv("data/weight_guesses.csv", index_col=0)
# drop some signatures since we have no associate mutationnal data
signatures = signatures.reindex(smoothed_mutation.index)

# which PCAs of smoothed mutations (mutated pathways?) 
# generate particular mutational signatures
loadings,_,_,_ = np.linalg.lstsq(transformed_features, signatures.values)

# 
import seaborn as sns
sns.clustermap(pd.DataFrame(loadings, columns=signatures.columns).T,cmap="seismic")

# replace the ensembl ids in the collumnns with hgnc gene symbols
smoothed_mutation_symbol = smoothed_mutation.copy()
smoothed_mutation_symbol.columns = smoothed_mutation.columns.map(lambda x : ensembl_symbol.get(x, "Pass"))

# 
sns.clustermap(smoothed_mutation_symbol, vmax=0.05, metric='cosine')
plt.savefig("figures/smoothing_"+parameter+"/clustermap_smoothed_mutation.png", dpi=500)
plt.savefig("figures/smoothing_"+parameter+"/clustermap_smoothed_mutation.svg")
