import numpy as np
import warnings

def createBins(adata, 
                bin_by, 
                min_cells, 
                target_totalBins = 1000, 
                method = ["mean", "median", "percentile"], 
                percentile = 75,
                use_raw = False):
    """
    Create pseudo-bulk bins for original data.
    """
    # get cell meta-data
    cell_meta = adata.obs

    # create bin groups
    cell_meta["bin_group"] = ""
    for i in bin_by:
        if i not in cell_meta.columns:
            raise Exception(i + " is not a valid cell annotation in anndata.")
        else:
            cell_meta["bin_group"] = [cell_meta["bin_group"][j] + "-" + cell_meta[i][j] for j in range(cell_meta.shape[0])]
    cell_meta["bin_group"] = [i.lstrip("-") for i in cell_meta["bin_group"]]

    # create binning factor
    cell_counts = cell_meta["bin_group"].value_counts()
    selected_bin_groups = cell_counts[cell_counts >= min_cells].index.values
    cell_meta = cell_meta[cell_meta["bin_group"].isin(selected_bin_groups)]

    cell_meta = cell_meta.sort_values("n_counts", ascending = False)
    cell_count_list = list(cell_counts)
    bin_factor = calculate_bin_size_factor(cell_numList = cell_count_list, target_totalBins = target_totalBins)

    # create bins using chunk2
    for bg in np.unique(cell_meta["bin_group"]):
        cells = cell_meta.loc[cell_meta["bin_group"] == bg, :].index.values
        if (len(cells) ** binfactor) < 2:
            n_bins = 2
        else:
            n_bins = round(len(cells) ** binfactor)
        cells_split = get_chunk_id(np.array_split(cells, n_bins))
        

        

def get_chunk_id(a):
    chunk_ids = []
    for i in range(len(a)):
        chunk_ids += [i] * len(a[i])
    return chunk_ids


def calculate_bin_size_factor(cell_numList, target_totalBins = 1000, bin_numVariation_allowed = 10):
    """
    Calculate binning size factor.

    Parameters
    ----------
    cell_numList
            A list of number of cells
    target_totalBins
            Final number of bins we want to get for the heatmap visualization
    bin_numVariation_allowed
            range of gap allowed between real and target number of bins
    """

    if len(cell_numList) > target_totalBins:
        raise Exception("Number of unique cell groups should be less than target total number of bins.")
    
    sum_cells = np.sum(cell_numList)
    if sum_cells < target_totalBins:
        warnings.warn("Number of total cells is less than target total number of bins.")

    Lower = 0.01; Upper = 1
    # perform Netibian bisection to converge upon an ideal nth-root
    for iteration in range(20):
        cf = (Lower + Upper) / 2
        sum = 0
        for n in cell_numList:
            # calculate number of bins per bin group (set min size = 2)
            if n ** cf < 2:
                num_bins = 2
            else:
                num_bins = round(n ** cf) # root value n ** cf is the number of bins in each cell label
        sum += num_bins
    
        # test if we need to raise or reduce the converging value
        if sum < (target_totalBins - bin_numVariation_allowed):
            Lower = cf
        elif sum > (target_totalBins + bin_numVariation_allowed):
            Upper = cf
        else:
            break

    return cf

