import random
from random import shuffle
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pylab import *

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

def create_cycler_colors(color_scheme):
    """
    create a list of colors extracted from a public
    """
    cmap = cm.get_cmap(color_scheme)    # PiYG
    cycler_colors = []

    for i in range(cmap.N):
        rgba = cmap(i)
        # rgb2hex accepts rgb or rgba
        cycler_colors.append(matplotlib.colors.rgb2hex(rgba)) 
    
    return cycler_colors

# create a color scheme using tab20b + tab20c
part_a = create_cycler_colors("tab20b")
part_b = create_cycler_colors("tab20c")
shuffle(part_a)
shuffle(part_b)

cycler_colors = part_a + part_b
plt.rcParams["axes.prop_cycle"] = cycler(color = cycler_colors)


def heatmap(shred, bin_type, save_output = True):
    """
    Draw the heatmap for gene modules
    """
    # format dataframes 
    df_subsetDEG = shred.shred_modules_df_2.copy() # load the gene module annotation dataframe

    if bin_type == "bin":
        df_heatmap = shred.heatmap_matrix.copy() # load heatmap matrix
        df_bin_meta = shred.bin_metadata.loc[df_heatmap.columns, :].copy() # order the bin metadata
    elif bin_type == "superbin":
        df_heatmap = shred.heatmap_matrix_super.copy()
        df_bin_meta = shred.superbin_metadata.loc[df_heatmap.columns, :].copy() # order the bin metadata
    else:
        raise Exception("bin_type should be either bin or superbin.")

    n_genes = df_heatmap.shape[0]
    n_bins = df_heatmap.shape[1]
    n_modules = len(np.unique(df_subsetDEG["Status"]))

    # initialize
    fig = plt.figure(figsize = (15, 15))
    gs = GridSpec(nrows = 3, ncols = 4, height_ratios=(0.10,0.05, 0.85), width_ratios = (0.70, 0.04, 0.07, 0.19))
    gs.update(wspace = 0.025, hspace = 0.025)

    #####----- bin metadata barplot ------#####
    ax0 = fig.add_subplot(gs[0,0])
    columns = []
    labels_for_legend = [] #initialize for legends
    colors_for_legend = []
    for column in df_bin_meta.columns:
        if column != "n_cells":
            columns.append(column)
            df_size = pd.DataFrame(df_bin_meta.groupby(columns).size())
            df_size.rename({0: "n_cells"}, axis = 1, inplace = True)
            labels = [i[-1] if type(i) == tuple else i for i in df_size.index.values]
            labels_for_legend += labels
            
            color_dict = dict(zip(labels, cycler_colors[:df_size.shape[0]])) # avoid different colors for duplicate keys
            count_aggregrated = 0
            for idx, label in enumerate(labels):
                count = df_size.iloc[idx, 0]
                ax0.barh(y = column, width = count, height = 1, left = count_aggregrated, color = color_dict[label])
                ax0.set_xticklabels("")
                ax0.set_xticks([])
                count_aggregrated += count

    # legend initialization
    colors = [matplotlib.colors.rgb2hex(ax0.patches[i].get_facecolor()) for i in range(len(ax0.patches))]
    handles = [plt.Rectangle((0,0),1,1, color = i) for i in colors]
    labels_unique = []
    handles_unique = []

    # remove duplicates
    for k in range(len(labels_for_legend)):
        if not labels_for_legend[k] in labels_unique:
            labels_unique.append(labels_for_legend[k])
            handles_unique.append(handles[k])

    ax0_1eg = fig.add_subplot(gs[0,3])
    leg1 = ax0_1eg.legend(handles_unique, labels_unique, prop = {'size': 10}, 
                          loc = 'upper left', title = "bin metadata", ncol = 5)
    leg1._legend_box.aligh = "left"
    ax0_1eg.add_artist(leg1)
    ax0_1eg.axis('off')

    # other formats
    ax0.set_xlim(0, n_bins)
    ax0.spines['left'].set_visible(False)
    ax0.spines['right'].set_visible(False) 
    ax0.spines['top'].set_visible(False) 
    ax0.spines['bottom'].set_visible(False) 

    #####----- n_cell histograms ------#####
    counts = df_bin_meta["n_cells"]
    ax1 = fig.add_subplot(gs[1,0])
    ax1.bar(x = list(range(n_bins)), height = counts, width = 1)
    ax1.set_xticklabels("")
    ax1.set_xticks([])
    ax1.set_xlim(0, n_bins)
    ax1.set_ylabel("n_cells", rotation = 0, labelpad = 24)

    #####----- gene module barplots ------#####
    ax2 = fig.add_subplot(gs[2,2])
    columns = ["Status", "reference_group","shred_plan"]
    labels_for_legend = []

    for column in columns:
        count_aggregrated = 0
        labels = np.unique(df_subsetDEG[column])
        labels_for_legend += list(labels)
        label_counts = dict(df_subsetDEG[column].value_counts())
        for label in labels:
            count = label_counts[label]
            ax2.bar(x = column, height = count, width = 1, bottom = count_aggregrated)
            ax2.set_yticklabels("")
            ax2.set_yticks([])
            count_aggregrated += count
    ax2.set_xticklabels(columns, rotation = 35, ha = "right")
    ax2.set_ylim(0, n_genes)
    ax2.spines['left'].set_visible(False) 
    ax2.spines['right'].set_visible(False) 
    ax2.spines['top'].set_visible(False) 
    ax2.spines['bottom'].set_visible(False) 

    # set legend
    n_Status = len(np.unique(df_subsetDEG["Status"]))
    n_reference_groups = len(np.unique(df_subsetDEG["reference_group"]))

    colors = [matplotlib.colors.rgb2hex(ax2.patches[i].get_facecolor()) for i in range(len(ax2.patches))]
    handles = [plt.Rectangle((0,0),1,1, color = i) for i in colors]
    # legend for status
    ax2_leg1 = fig.add_subplot(gs[2,3])
    leg1 = ax2_leg1.legend(handles[:n_Status], labels_for_legend[:n_Status], 
                            prop = {'size': 10}, loc = 'upper left', title = "Status", ncol=4)
    leg1._legend_box.aligh = "left"
    ax2_leg1.add_artist(leg1)
    # legend for others (reference and shred_plan)
    ax2_leg2 = fig.add_subplot(gs[2,3])
    leg2 = ax2_leg2.legend(handles[n_Status:(n_Status + n_reference_groups)], labels_for_legend[n_Status:(n_Status + n_reference_groups)], 
                            prop = {'size': 10}, loc = 'center left',title = "reference group")
    leg3 = ax2_leg2.legend(handles[(n_Status + n_reference_groups):], labels_for_legend[(n_Status + n_reference_groups):], 
                            prop = {'size': 10}, loc = 'lower left',title = "shred plan")
    leg2._legend_box.aligh = "left"
    leg3._legend_box.aligh = "left"
    ax2_leg2.add_artist(leg2)
    ax2_leg2.add_artist(leg3)

    ax2_leg1.axis('off')
    ax2_leg2.axis('off')

    #####----- gene module score histogram ------##### 
    ax3 = fig.add_subplot(gs[2,1])
    ax3.barh(y = list(range(n_genes)), width = list(df_subsetDEG["logfoldchanges"]), height = 1)
    ax3.xaxis.tick_top()
    ax3.set_yticklabels("")
    ax3.set_yticks([])
    ax3.set_ylim(0, n_genes)
    ax3.set_xlabel("logFC", rotation = 35, ha = "right")
    ax3.plot([1, 1], [0, 5200], color = "#990000", linestyle = "dashed", lw = 0.8)

    #####----- main heatmap ------#####
    ax5 = fig.add_subplot(gs[2,0])
    sns.heatmap(df_heatmap, vmin = 0, vmax = 10, yticklabels=False, xticklabels=False, cmap = "bwr", ax = ax5, cbar = False)
    ax5.set_ylabel("")

    output_name = (shred.output_folder + "/figures/heatmap.png") if bin_type == "bin" else (shred.output_folder + "/figures/heatmap_superbin.png")
    if save_output:
        fig.savefig(output_name, bbox_inches = "tight")


def correlation_bins(shred):
    """
    Draw the correlation map for bins.
    """
    return 0


