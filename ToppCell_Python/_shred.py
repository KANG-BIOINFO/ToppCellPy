import numpy as np
import pandas as pd 
import scanpy as sc
from .differential import compute_differential_analysis

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
        order_bins,
        order_modules = None,
        method = "wilcoxon"
    ):
        __init__()
        self.adata = adata
        self.shred_plan = shred_plan
        self.bin_group = bin_group
        self.order_bins = order_bins
        self.order_modules = order_modules
        self.method = method
    
    def do_shredplan(self):
        """
        Run the user-customed shred plan
        """
        for sub_plan in self.shred_plan:
            # get target and reference
            if "|" not in sub_plan:
                target = sub_plan; reference = None
            else:
                target = sub_plan.split("|")[0]; reference = sub_plan.split("|")[1]
            df_deg = compute_levelWise_differential_analysis(self, target = target, reference = reference)

            self.shred_module[sub_plan] = df_deg
    
    def create_heatmap_matrix(self):
        """
        Create heatmap matrix.
        """
        


