import pandas as pd 
import numpy as np 
import scanpy as sc
import pickle
import ToppCellPy as tp
import time 

path = "/Users/kang/Dropbox/Code/data/toppcell_test/"

adata = sc.read("/Users/kang/Dropbox/Code/data/batch2_all_normalized_filtered_subset.h5ad")
print(adata)
# set up shred object and plan
shred = tp.Shred(adata = adata,
            shred_plan = ["stim", "cell", "stim+cell|stim"],
            bin_group = ["stim", "cell"],
            order_bins = None,
            order_modules = None,
            method = "wilcoxon",
            output_dir = "/Users/kang/Dropbox/Code/data/toppcell_test/")

# run everything once.
# shred.toppcell_batchRun(enrich_categories = ["GeneOntologyCellularComponent"])

# run levelwise DE analysis and generate gene modules
shred.do_shredplan()

# make heatmap table for visualization
shred.create_heatmap_matrix()

# draw heatmap figure
shred.draw_heatmap()

# do enrichment for all modules
shred.enrich_modules(categories = ["GeneOntologyCellularComponent"])

# draw toppcluster plot
shred.toppcluster()

# create GCT files
shred.createGCT()