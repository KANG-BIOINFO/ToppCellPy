import pandas as pd
import numpy as np
import gseapy as gp 
import pickle
from os import path
import pkg_resources

THETA_LOWER_BORDER = 1e-200
DATA_PATH = pkg_resources.resource_filename('ToppCellPy', '').rstrip("ToppCellPy") + "data/ToppGene_ref/"

def module_enrich_ranked(ranked_gene_table, 
                    terms = ["GeneOntologyMolecularFunction",
                                "GeneOntologyBiologicalProcess",
                                "GeneOntologyCellularComponent",
                                "MousePheno",
                                "Pathway"],
                    min_size = 1, max_size = 1500):
    """
    Enrichment for ranked gene list.
    """

    # load all reference gene sets (dictionaries)
    enrich_all_dicts = {}
    for term in terms:
        with open (DATA_PATH + term + ".pl", "rb") as f:
            enrich_dict = pickle.load(f)
        enrich_all_dicts[term] = enrich_dict

    # do enrichment for each category
    gsea_output = pd.DataFrame()
    for term in terms:
        ref = enrich_all_dicts[term]
        pre_res = gp.prerank(rnk = ranked_gene_table,
                            gene_sets = ref,
                            processes = 4, permutation_num = 100,
                            min_size = min_size, max_size = max_size, outdir = None, no_plot = True)
        res2d = pre_res.res2d .sort_values(["fdr"], ascending = True)
        res2d["Category"] = term
        gsea_output = pd.concat([gsea_output, res2d], axis = 0)

    return gsea_output 

def module_enrich(gene_table, 
                    terms = ["GeneOntologyMolecularFunction",
                                "GeneOntologyBiologicalProcess",
                                "GeneOntologyCellularComponent",
                                "MousePheno",
                                "Pathway"]):
    """
    Enrichment for gene list.
    """

    # load all reference gene sets (dictionaries)
    enrich_all_dicts = {}
    for term in terms:
        with open (DATA_PATH + term + ".pl", "rb") as f:
            enrich_dict = pickle.load(f)
        enrich_all_dicts[term] = enrich_dict

    # do enrichment for each category
    gsea_output = pd.DataFrame()
    for term in terms:
        ref = enrich_all_dicts[term]
        enrichr_res = gp.enrichr(gene_list = gene_table,
                                gene_sets = ref,
                                cutoff = 1,
                                outdir = None, no_plot = True)
        res2d = enrichr_res.res2d .sort_values(["Adjusted P-value"], ascending = True)
        res2d["Category"] = term
        gsea_output = pd.concat([gsea_output, res2d], axis = 0)

    return gsea_output     

def apply_toppcluster(shred):
    """
    Apply ToppCluster enrichment across gene modules.
    """
    df_module_enrichment = shred.df_module_enrichment.copy()
    df_module_enrichment["enrichment score"] = [-np.log10(i + THETA_LOWER_BORDER) for i in df_module_enrichment["Adjusted P-value"]]
    
    all_enriched_terms = np.unique(df_module_enrichment["Term"])
    all_modules = np.unique(df_module_enrichment["module"])

    df_toppcluster_neglog10pval = pd.DataFrame(index = all_modules, columns = all_enriched_terms, data = 0)
    
    for row in range(df_module_enrichment.shape[0]):
        term = df_module_enrichment.iloc[row, 1]
        module = df_module_enrichment.iloc[row, 7]
        score = df_module_enrichment.iloc[row, 8]
        df_toppcluster_neglog10pval.loc[module, term] = score 

    return df_toppcluster_neglog10pval
    
def test_load_data():
    for term in ["GeneOntologyMolecularFunction"]:
        with open (DATA_PATH + term + ".pl", "rb") as f:
            enrich_dict = pickle.load(f)
            print(enrich_dict)


    
