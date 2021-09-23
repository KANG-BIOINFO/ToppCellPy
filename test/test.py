import pandas as pd 
import numpy as np 
import scanpy as sc
import sys

sys.path.insert(1, '/Users/kang/Dropbox/Code/ToppCell-Python/')
import ToppCell_Python as tp

path = "/Users/kang/Dropbox/Code/data/toppcell_test/"

adata = sc.read("/Users/kang/Dropbox/Code/data/batch2_all_normalized_filtered.h5ad")

# create bins
'''
df_bin_meta, df_bins = tp.createBins(adata, bin_by = ["stim", "cell"], min_cells = 5, target_totalBins = 200)
df_bins.to_csv(path + "test_bin_table.txt", sep = "\t")
df_bin_meta.to_csv(path + "test_bin_meta_table.txt", sep = "\t")
'''

# test shred class
shred = tp.Shred(adata = adata,
            shred_plan = ["stim", "cell", "stim+cell|stim"],
            bin_group = ["stim", "cell"],
            order_bins = None,
            order_modules = None,
            method = "wilcoxon")
df_deg_combined = shred.do_shredplan()
df_deg_combined.to_csv(path + "deg_shredplan.txt", sep = "\t")
df_heatmap = shred.create_heatmap_matrix(path + "heatmap_matrix", sep = "\t")

# make heatmap table


# generate heatmap figure 