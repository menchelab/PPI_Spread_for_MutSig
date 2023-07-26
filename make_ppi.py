import pandas as pd
import networkx as nx
import numpy as np

# load the edge list
df = pd.read_csv("data/9606.protein.physical.links.v11.5.txt", sep=" ")

# set a cutoff for score
cutoff_score = 600
df = df[df.combined_score > cutoff_score]

ensembl_map = pd.read_csv("data/Homo_sapiens.GRCh38.110.entrez.tsv", sep="\t").set_index("protein_stable_id")["gene_stable_id"].to_dict()
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

nx.write_gml(ppi, "data/ppi.gml")


pos = nx.spring_layout(ppi)
pd.DataFrame(pos).T.to_csv("data/ppi_layout.csv")
