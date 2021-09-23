import pandas as pd 
import numpy as np 
import scanpy as sc
import sys

sys.path.insert(1, '/Users/kang/Dropbox/Code/ToppCell-Python/')
import ToppCell_Python as tp

path = "/Users/kang/Dropbox/Code/data/toppcell_test/"

adata = sc.read("/Users/kang/Dropbox/Code/data/batch2_all_normalized_filtered.h5ad")

# create bins
df_bins = tp.createBins(adata, bin_by = ["stim", "cell"], min_cells = 5, target_totalBins = 200)
df_bins.to_csv(path + "test_bin_table.txt", sep = "\t")

# 