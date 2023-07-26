
import pandas as pd
import json
import sklearn.decomposition as sd
import numpy as np
import seaborn as sns

import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt

import os
import sys
parameter =  sys.argv[1]

transformed_features = pd.read("data/pca_smooth_features_"+parameter+".csv")
# helper dictionary to map ensembl to hgnc symbol
name_mapper = json.load(open("data/protein-coding_gene.json"))['response']['docs']ensembl_symbol = {i['ensembl_gene_id']:i['symbol'] for i in name_mapper if 'ensembl_gene_id' in i.keys()}


import umap
um = umap.UMAP(n_components=2)
xy = um.fit_transform(transformed_features)
fig,ax = plt.subplots()
ax.scatter(xy[:,0], xy[:,1])

fig.savefig("figures/smoothing_"+parameter+"/umap.svg")
fig.savefig("figures/smoothing_"+parameter+"/umap.png", dpi=500)

# 
os.makedirs("figures/smoothing_"+parameter+"/pca/", exist_ok = True) 
i = 0
j = 1
fig,ax = plt.subplots()
ax.scatter(transformed_features.loc[:,i], transformed_features.loc[:,j])
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

