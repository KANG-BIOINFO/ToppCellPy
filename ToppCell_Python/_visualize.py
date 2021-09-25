import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import random
from scipy.stats import pearsonr

params = {
        "image.cmap": "tab20c",
        "font.size": 16,
        "font.family": "Arial",
        "axes.titlesize": 18,
        "axes.labelsize": 14,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
        "legend.title_fontsize": 16,
        "figure.figsize": (8, 8),
        "axes.xmargin": 0.1,
        "axes.ymargin": 0.04}
plt.rcParams.update(params)


def heatmap(shred, output_name = "heatmap.png"):
    """
    Draw the heatmap for gene modules
    """
    # format dataframes
    df_heatmap = shred.heatmap_matrix.copy() # load heatmap matrix
    df_bin_meta = shred.bin_metadata.loc[df_heatmap.columns, :].copy() # order the bin metadata
    df_subsetDEG = shred.shred_modules_df_2.copy() # load the gene module annotation dataframe
    n_genes = df_heatmap.shape[0]
    n_bins = df_heatmap.shape[1]
    n_modules = len(np.unique(df_subsetDEG["Status"]))

    # initialize
    fig = plt.figure(figsize = (10,10))
    gs = GridSpec(nrows = 3, ncols = 3, height_ratios=(0.10,0.05, 0.85), width_ratios = (0.85,0.05,0.10))
    gs.update(wspace = 0.025, hspace = 0.025)

    # bin metadata barplot
    ax3 = fig.add_subplot(gs[0,0])
    columns = df_bin_meta.columns

    for column in columns:
        count_aggregrated = 0
        labels = np.unique(df_bin_meta[column])
        label_counts = dict(df_bin_meta[column].value_counts())
        for label in labels:
            count = label_counts[label]
            ax3.barh(y = column, width = count, height = 1, left = count_aggregrated)
            ax3.set_xticklabels("")
            ax3.set_xticks([])
            count_aggregrated += count
    # ax3.set_yticklabels(columns, rotation = 45)
    ax3.set_xlim(0, n_bins)
    ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_visible(False) 
    ax3.spines['top'].set_visible(False) 
    ax3.spines['bottom'].set_visible(False) 

    # n_cell histograms
    counts = np.random.normal(100, 20, n_bins) # need to get back and fix
    ax0 = fig.add_subplot(gs[1,0])
    ax0.bar(x = list(range(n_bins)), height = counts, width = 1)
    ax0.set_xticklabels("")
    ax0.set_xticks([])
    ax0.set_xlim(0, n_bins)
    ax0.set_ylabel("n_cells")

    # gene module barplots
    ax1 = fig.add_subplot(gs[2,2])
    columns = ["Status", "reference_group","shred_plan"]

    for column in columns:
        count_aggregrated = 0
        labels = np.unique(df_subsetDEG[column])
        label_counts = dict(df_subsetDEG[column].value_counts())
        for label in labels:
            count = label_counts[label]
            ax1.bar(x = column, height = count, width = 1, bottom = count_aggregrated)
            ax1.set_yticklabels("")
            ax1.set_yticks([])
            count_aggregrated += count
    ax1.set_xticklabels(columns, rotation = 45)
    ax1.set_ylim(0, n_genes)
    ax1.spines['left'].set_visible(False) 
    ax1.spines['right'].set_visible(False) 
    ax1.spines['top'].set_visible(False) 
    ax1.spines['bottom'].set_visible(False) 

    # gene module score histogram
    ax4 = fig.add_subplot(gs[2,1])
    ax4.barh(y = list(range(n_genes)), width = list(df_subsetDEG["logfoldchanges"]), height = 1)
    ax4.set_yticklabels("")
    ax4.set_yticks([])
    ax4.set_ylim(0, n_genes)
    ax4.set_xlabel("logFC")

    # heatmap itself
    ax2 = fig.add_subplot(gs[2,0])
    sns.heatmap(df_heatmap, vmin = 0, vmax = 10, yticklabels=False, xticklabels=False, cmap = "bwr", ax = ax2, cbar = False)

    fig.savefig(output_name)


def draw_module_enrichment(shred, top_n_modules):
    df_module_enrichment = shred.df_module_enrichment


def correlation_bins(shred):
    """
    Draw the correlation map for bins.
    """
    return 0