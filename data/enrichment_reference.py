import pandas as pd 
import numpy as np 
import pickle 
import os 

path = "/Users/jinmr2/Dropbox/Code/data/ToppGene_Data/"
files = [i for i in os.listdir(path) if i.endswith(".txt")]
for enrich_file in files:
    name = enrich_file.rstrip(".txt")
    df = pd.read_csv(path + enrich_file, sep = "\t", header = 0, index_col = 0)
    enrich_dict = {}
    for gene_set in np.unique(df["concept_name"]):
        genes = df.loc[df["concept_name"] == gene_set, "symbol"]
        enrich_dict[gene_set] = list(genes)
    with open(path + name + ".pl", "wb") as f:
        pickle.dump(enrich_dict, f)
