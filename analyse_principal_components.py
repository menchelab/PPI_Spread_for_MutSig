import networkx as nx
import pandas as pd
import sys

import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt

parameter = sys.argv[1]
pcs = pd.read_csv("data/pca_components_"+parameter+".csv", index_col=0)

layout = pd.read_csv("data/ppi_layout.csv", index_col=0)
ppi = nx.read_gml("data/ppi.gml")

import seaborn as sns


aaa = sns.clustermap(pcs.T, cmap="seismic", metric = 'cosine', col_cluster=False, center=0, vmax = 0.1, vmin=-0.1)
aaa.savefig("figures/smoothed_"+parameter+"/pc_clustermap,png")
aaa.fig.show()

for pc in pcs.columns[:10]:
    fig,ax = plt.subplots(figsize=[10,10])
    nx.draw_networkx_nodes(ppi, ax=ax, pos = {i:v for i,v in layout.iterrows()}, node_color=pcs[pc].as_dict())
    nx.draw_networkx_edges(ppi, ax=ax, pos = {i:v for i,v in layout.iterrows()})
    fig.savefig("figures/smoothed_"+parameter+"/pc.png", dpi=500)


