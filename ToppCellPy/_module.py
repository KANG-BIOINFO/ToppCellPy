import numpy as np
import pandas as pd 
import scanpy as sc

class ToppCellModule:
    """
    ToppCell container for each gene module. (recommand for all genes)

    Parameters
    ----------
    label_level
        Cell label (annotation level) for the comparison.
    target_label
        Cell label for target cells.
    reference_label
        Cell label for reference cells.
    method
        Statistical method used for the significance test.
    genes
        List for genes for the module.
    pvals
        Adjusted p values for the significance of genes in the list.
    logfoldchanges
        Log fold changes for the gene list.
    pts
        Percentage of cells with non-zero expression level within group for each gene.
    pts_rest
        Percentage of cells with non-zero expression level out of group for each gene.
    """
    def __init__(
        self,
        label_level,
        target_label,
        reference_label,
        method,
        genes,
        pvals,
        logfoldchanges,
        pts,
        pts_rest,
    ):
        __init__()
        self.label_level = label_level 
        self.target_label = target_label
        self.reference_label = reference_label
        self.method = method
        self.genes = genes
        self.pvals = pvals
        self.logfoldchanges = logfoldchanges
        self.pts = pts
        self.pts = pts_rest

        self.n_size = len(self.genes)
        
    def enrichment(self, top_n_genes, categories, pval_adjust_method, pval_cutoff, min_genes_pathway, max_genes_pathway):
        """
        Do the enrichment for the gene list.

        Parameters
        ----------
        categories
            Categories of enrichment analysis
        pval_adjust_method
            P value adjustment method
        pval_cutoff
            P value cutoff to get significant enrichment
        min_genes_pathway
            minimal number of gens in reference pathway
        max_genes_pathway
            maximal number of genes in reference pathway
        """
        
    def cell_type_projection(self, )

