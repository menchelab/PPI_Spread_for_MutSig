import networkx as nx
import pandas as pd
import numpy as np
import sys

parameter = sys.argv[1]
pcs = pd.read_csv("data/pca_smooth_features_"+parameter+".csv", index_col=0)

import gseapy as gp

pcs.columns = pcs.columns.astype(int)

out = []
for pc in pcs.columns[:5]:
    gene_list = pcs.index[pcs[pc] > np.quantile(pcs[pc], 0.95)]

    import json
    name_mapper = json.load(open("data/protein-coding_gene.json"))['response']['docs']
    ensembl_symbol = {i['ensembl_gene_id']:i['symbol'] for i in name_mapper if 'ensembl_gene_id' in i.keys()}

    full_gene_list = list(set(pcs.index) & set(ensembl_symbol.keys()))
    full_gene_list = [ensembl_symbol[i] for i in full_gene_list]

    gene_list = list(set(gene_list) & set(ensembl_symbol.keys()))
    gene_list = [ensembl_symbol[i] for i in gene_list]

    enr = gp.enrichr(gene_list=gene_list,
                      gene_sets=['GO_Biological_Process_2023'],
                      organism='human', 
                      outdir=None,
                      background=full_gene_list
                     )

    df = enr.results
    df["log_p"] = df["Adjusted P-value"].map(lambda x:-np.log(x))
    df["GO Name"] = df["Term"].map(lambda x: x.split("(")[0].strip())
    df["PC"] = pc
    out.append(df)

go_term_associations = pd.concat(out)

