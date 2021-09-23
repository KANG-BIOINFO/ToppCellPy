import pandas as pd 
import numpy as np 
import scanpy as sc
import sys
sys.path.insert(1, '/Users/jinmr2/Dropbox/Code/ToppCell-Python/')

import ToppCell_Python as tp
adata = sc.read("/Users/jinmr2/Dropbox/Code/data/batch2_all_raw_filtered.h5ad")
df_bins = tp.createBins(adata, bin_by = ["stim", "cell"], min_cells = 5, target_totalBins = 200)
print(df_bins)
df_bins.to_csv("test_bin_table.txt", sep = "\t")