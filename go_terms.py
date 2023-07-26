import pandas as pd
import json

# Read the gene annotations from the file
go_annotations = pd.read_csv("goa_human.gaf", skiprows=41, sep="\t", header=None)
# Specify the column indices for gene and GO terms
gene_column = 2
go_column = 4



# Load the GO ontology file
go_info = json.load(open("go-basic.json"))
# Create a dictionary mapping GO IDs to their corresponding names
go_id_to_name = {i["id"].split("/")[-1].replace("_", ":"): i['lbl'] for i in go_info['graphs'][0]['nodes'] if 'lbl' in i.keys()}


# Create a dictionary to store gene annotations
annotations = {}
# Process each gene column in the dataframe
for gene in df.columns:
    # Select rows where the gene column matches the current gene
    # NB: if you have different gene identifiers in your own dataset
    # you will have to adjust this line to select the right GO terms
    selected_rows = go_annotations[gene_column] == gene

    # Get unique GO IDs associated with the current gene
    go_ids = set(go_annotations[selected_rows][go_column])

    # Map GO IDs to their corresponding names and store them in the annotations dictionary
    annotations[gene] = [go_id_to_name[go_id] for go_id in go_ids]