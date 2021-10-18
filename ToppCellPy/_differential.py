import numpy as np
import pandas as pd
import scanpy as sc

def compute_levelWise_differential_analysis(shred, target, reference, plan_name, module_name_delimiter = "~"):
    """
    Do differential analysis.

    Parameters
    ----------
    shred
            Shred object
    target
            numerator level for DE analysis, e.g. stim or stim+cell
    reference
            denominator level for DE analysis, e.g. None (World)
    plan_name
            name of shred plan
    module_name_delimiter
            delimiter character in module names. e.g. "COVID-19~Monocyte"
    """
    adata = shred.adata.copy()
    cell_meta = adata.obs

    # create a new column for target and reference
    if "+" in target:
        levels_target = target.split("+")
        for idx, level in enumerate(levels_target):
            if idx == 0:
                target_value = cell_meta[level]
            else:
                target_value = [(target_value[i] + "~" + cell_meta[level][i]) for i in range(cell_meta.shape[0])]
    else:
        levels_target = target
        target_value = cell_meta[target]
    cell_meta["target_value"] = target_value

    if reference != None:
        if "+" in reference:
            levels_ref = reference.split("+")
            for idx, level in enumerate(levels_ref):
                if idx == 0:
                    reference_value = cell_meta[level]
                else:
                    reference_value = [(reference_value[i] + "-" + cell_meta[level][i]) for i in range(cell_meta.shape[0])]
        else:
            reference_value = cell_meta[reference]
        cell_meta["reference_value"] = reference_value
    adata.obs = cell_meta

    # compute DEG
    # when the reference is global (omitted)
    if reference == None:
        sc.tl.rank_genes_groups(adata, groupby = "target_value", pts = True, method = shred.method)
        df_deg = format_DEGs(adata)
        df_deg["reference_group"] = ""
        df_deg["shred_plan"] = plan_name
    else:
        df_deg = pd.DataFrame()
        for reference_i in np.unique(adata.obs["reference_value"]):
            adata_sub = adata[adata.obs["reference_value"] == reference_i, : ].copy()
            sc.tl.rank_genes_groups(adata_sub, groupby = "target_value", pts = True, method = shred.method)
            df_deg_sub = format_DEGs(adata_sub)
            df_deg_sub["reference_group"] = reference_i
            df_deg = pd.concat([df_deg, df_deg_sub], axis = 0)
        df_deg["shred_plan"] = plan_name

    # add target level metadata for each module (decompose status value)
    for idx, level in enumerate(levels_target):
        df_deg[level] = [i.split(module_name_delimiter)[idx] for i in df_deg["Status"]]
        
    return df_deg


def format_DEGs(adata):
    """
    Format the output from anndata
    """
    keys = ["names","scores","logfoldchanges","pvals","pvals_adj","pts","pts_rest"]
    for i,key in enumerate(keys):
        a = pd.DataFrame(adata.uns["rank_genes_groups"][key]) # transfer to data frame
        b = pd.DataFrame(a.values.T.reshape(1,a.shape[0]*a.shape[1]).T) # reformat the data frame to one column
           
        if i == 0:
            b.columns = [key] # rename the column name
            b["Status"] = sorted(list(a.columns)*a.shape[0]) # add Status annotation
            b.set_index([key],inplace=True)
            b_merged = b
        else:
            if key in ["pts","pts_rest"]:
                pts_all = []
                for cell_group in np.unique(b_merged["Status"]):
                    genes = b_merged.loc[b_merged["Status"] == cell_group,:].index.values
                    pts_all = pts_all + list(a.loc[genes, cell_group])
                b_merged[key] = pts_all
            else:
                b_merged[key] = list(b[0])
        
    return b_merged