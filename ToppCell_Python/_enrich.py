import pandas as pd
import numpy as np
import gseapy as gp 
import pickle
from os import path

def module_enrich(ranked_gene_table, 
                  terms = ["GeneOntologyMolecularFunction",
                            "GeneOntologyBiologicalProcess",
                            "GeneOntologyCellularComponent",
                            "MousePheno",
                            "Pathway"],
                  min_size = 5, max_size = 1500):
    # load all reference gene sets (dictionaries)
    loc = "/Users/jinmr2/Dropbox/Code/data/ToppGene_Data/" if path.isdir("/Users/jinmr2/Dropbox/Code/data/ToppGene_Data/") else "/Users/kang/Dropbox/Code/data/ToppGene_Data/"
    enrich_all_dicts = {}
    for term in terms:
        with open (loc + term + ".pl", "rb") as f:
            enrich_dict = pickle.load(f)
        enrich_all_dicts[term] = enrich_dict

    # do enrichment for each category
    gsea_output = pd.DataFrame()
    for term in terms:
        ref = enrich_all_dicts[term]
        pre_res = gp.prerank(rnk = ranked_gene_table,
                            gene_sets = ref,
                            processes = 4, permutation_num = 100,
                            min_size = min_size, max_size = max_size)
        res2d = pre_res.res2d .sort_values(["fdr"], ascending = True)
        res2d["Category"] = term
        gsea_output = pd.concat([gsea_output, res2d], axis = 0)

    return gsea_output 

    
