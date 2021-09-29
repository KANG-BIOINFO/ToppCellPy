import pandas as pd 
import numpy as np 
import scanpy as sc
import sys
import pickle

sys.path.insert(1, '/Users/jinmr2/Dropbox/Code/ToppCell-Python/')
import ToppCell_Python as tp
import time 

path = "/Users/jinmr2/Dropbox/Code/data/toppcell_test/"

adata = sc.read("/Users/jinmr2/Dropbox/Code/data/batch2_all_normalized_filtered.h5ad")

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

# make heatmap table
df_heatmap, df_genemodule = shred.create_heatmap_matrix()
df_heatmap.to_csv(path + "heatmap_matrix.txt", sep = "\t")
df_genemodule.to_csv(path + "heatmap_geneModules.txt", sep = "\t")
shred.bin_metadata.to_csv(path + "bin_metadata.txt",sep = "\t")

# draw heatmap figure
shred.draw_heatmap(output_name = path + "heatmap.png")

# do enrichment for all modules
df_module_enrich = shred.enrich_modules(categories = ["GeneOntologyCellularComponent"])
df_module_enrich.to_csv(path + "module_enrichments.txt", sep = "\t")

# draw toppcluster plot
df_toppcluster = shred.toppcluster()
df_toppcluster.to_csv(path + "toppcluster_modules.txt", sep = "\t")