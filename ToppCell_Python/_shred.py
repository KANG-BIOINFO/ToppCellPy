import numpy as np
import pandas as pd 
import scanpy as sc
from ._differential import compute_levelWise_differential_analysis
from ._pseudo import createBins

class Shred:
    """
    An object of a shred job, including anndata, bin table, gene module lists and shred plan.

    Parameters
    ----------
    adata
            anndata for shred
    shred_plan
            shred plan for anndata. e.g. ["stim", "cell", "stim+cell|stim"]
    bin_group
            groups for binning plan in visualization
    order_bins
            orders of bins
    order_modules
            orders of modules
    method
            statistical methods for differential expression analysis
    """
    def __init__(
        self,
        adata,
        shred_plan,
        bin_group,
        bin_num = 1000,
        bin_min_cells = 5,
        order_bins = None,
        order_modules = None,
        method = "wilcoxon"
    ):
        self.adata = adata
        self.shred_plan = shred_plan
        self.bin_group = bin_group
        self.bin_num = 1000
        self.bin_min_cells = 5
        self.order_modules = order_modules
        self.method = method
        if order_bins == None:
            self.order_bins = self.bin_group
        else:
            self.order_bins = order_bins

        self.shred_module = {}
        self.module_groups = get_all_terms(shred_plan)
        
        self.bin_metadata, self.bin_matrix = createBins(adata, bin_by = bin_group, min_cells = bin_min_cells, target_totalBins = bin_num)
    
    def do_shredplan(self):
        """
        Run the user-customed shred plan
        """
        df_deg_combined = pd.DataFrame()
        for sub_plan in self.shred_plan:
            # get target and reference
            if "|" not in sub_plan:
                target = sub_plan; reference = None
            else:
                target = sub_plan.split("|")[0]; reference = sub_plan.split("|")[1]
            df_deg = compute_levelWise_differential_analysis(self, target = target, reference = reference, plan_name = sub_plan)

            self.shred_module[sub_plan] = df_deg
            df_deg_combined = pd.concat([df_deg_combined, df_deg], axis = 0)            
        
        return df_deg_combined
    
    def create_heatmap_matrix(self, plans, top_n_genes = 200):
        """
        Create heatmap matrix.
        """
        df_bin_meta = self.bin_metadata
        df_bin = self.bin_matrix
        
        # get bin table
        df_bin_meta = df_bin_meta.sort_values(self.bin_group)
        df_bin = df_bin.loc[list(df_bin_meta.columns), : ]

        # get organized gene modules ready
        df_DEG = compute_differential_analysis(self)
        df_subsetDEG = df_DEG.groupby(["Status"]).head(top_n_genes)
        df_subsetDEG = df_subsetDEG.sort_values(["shred_plan", "reference_group", "Status", "pts"], ascending = [True, True, True, False])

        # create heatmap
        df_heatmap = df_bin.loc[list(df_subsetDEG.index.values),:]
        self.heatmap_matrix = heatmap

        return heatmap

    def createJson(self):
        """
        create json file
        """
        #self.json = json
        return 1
    
    def createGCT(self):
        """
        create GCT file
        """
        
        return 1


def get_all_terms(terms):
    all_terms = []
    for term in terms:
        if "+" not in term:
            all_terms.append(term)
        else:
            all_terms += term.split("+")
            




        
        


