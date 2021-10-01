import os
import datetime
import numpy as np
import pandas as pd 
import scanpy as sc

import matplotlib.pyplot as plt 
import seaborn as sns

from ._differential import compute_levelWise_differential_analysis
from ._pseudo import createBins, createSuperbins
from ._visualize import heatmap
from ._enrich import module_enrich_ranked, module_enrich, apply_toppcluster

class Shred:
    """
    An object of a shred job, including anndata, bin table, gene module lists and shred plan.

    Parameters
    ----------
    adata
            anndata for shred
    shred_plan
            shred plan for anndata. e.g. ["stim", "cell", "stim+cell|stim"]
    bin_group
            groups for binning plan in visualization
    order_bins
            orders of bins
    order_modules
            orders of modules
    method
            statistical methods for differential expression analysis
    """

    def __init__(
        self,
        adata,
        shred_plan,
        bin_group,
        bin_num = 1000,
        bin_min_cells = 5,
        order_bins = None,
        order_modules = None,
        method = "wilcoxon",
        output_dir = "./",
        save_output = True
    ):
        self.adata = adata
        self.shred_plan = shred_plan # plan for shred (gene module generation)
        self.bin_group = bin_group # the way to make pseudo-bulk bins
        self.bin_num = 1000
        self.bin_min_cells = 5
        self.order_modules = order_modules # the order of modules in heatmap
        self.method = method
        self.save_output = save_output
        
        self.order_bins = self.bin_group if order_bins == None else order_bins
        self.shred_module = {}
        self.module_groups = get_all_terms(shred_plan)
        
        # create bins for heatmap visualization
        self.bin_metadata, self.bin_matrix = createBins(adata, bin_by = bin_group, 
                                                        min_cells = bin_min_cells, target_totalBins = bin_num)
        self.superbin_metadata, self.superbin_matrix = createSuperbins(adata, bin_by = bin_group)
        
        if output_dir != None:
            if not os.path.isdir(output_dir + "/output"):    
                self.output_folder = output_dir + "/output/"
            else:
                self.output_folder = output_dir + "/output_" + str(datetime.datetime.now()) + "/"
            os.mkdir(self.output_folder)
            
            if self.save_output:
                self.bin_metadata.to_csv(self.output_folder + "bin_metadata.txt", sep = "\t")
                self.bin_matrix.to_csv(self.output_folder + "bin_matrix.txt", sep = "\t")
                self.superbin_metadata.to_csv(self.output_folder + "superbin_metadata.txt", sep = "\t")
                self.superbin_matrix.to_csv(self.output_folder + "superbin_matrix.txt", sep = "\t")


    def do_shredplan(self):
        """
        Run the user-customized shred plan, where hierarchical gene modules will be generated.
        """
        df_deg_combined = pd.DataFrame()
        for sub_plan in self.shred_plan:
            # get target and reference
            if "|" not in sub_plan:
                target = sub_plan; reference = None
            else:
                target = sub_plan.split("|")[0]; reference = sub_plan.split("|")[1]
            df_deg = compute_levelWise_differential_analysis(self, target = target, reference = reference, plan_name = sub_plan)

            self.shred_module[sub_plan] = df_deg
            df_deg_combined = pd.concat([df_deg_combined, df_deg], axis = 0)            
        
        self.shred_modules_df = df_deg_combined        
        if self.save_output:
            self.shred_modules_df.to_csv(self.output_folder + "genemodules_all.txt", sep = "\t")
                
        # return df_deg_combined
    

    def create_heatmap_matrix(self, top_n_genes = 200):
        """
        Create heatmap matrix based on the generated gene modules.

        Parameters
        ----------
        top_n_genes
                Number of most significant genes to show in each gene module.
        """

        self.top_n_genes = top_n_genes

        df_bin_meta = self.bin_metadata
        df_bin = self.bin_matrix
        df_superbin_meta = self.superbin_metadata
        df_superbin = self.superbin_matrix
        
        # get bin table
        df_bin_meta = df_bin_meta.sort_values(self.bin_group)
        df_bin = df_bin[list(df_bin_meta.index.values)]
        df_superbin_meta = df_superbin_meta.sort_values(self.bin_group)
        df_superbin = df_superbin[list(df_superbin_meta.index.values)]

        # get organized gene modules ready
        try:
            df_DEG = self.shred_modules_df
        except AttributeError:
            raise Exception("Level-wise comparison should be done first using do_shredplan prior to creating heatmap.")

        df_subsetDEG = df_DEG.groupby(["Status"]).head(top_n_genes)
        df_subsetDEG = df_subsetDEG.sort_values(["shred_plan", "reference_group", "Status", "pts"], ascending = [True, True, True, False])

        # create heatmap
        df_heatmap = df_bin.loc[list(df_subsetDEG.index.values),:]
        df_heatmap = df_heatmap.astype(float)
        
        df_heatmap_super = df_superbin.loc[list(df_subsetDEG.index.values),:]
        df_heatmap_super = df_heatmap_super.astype(float)

        self.heatmap_matrix = df_heatmap
        self.heatmap_matrix_super = df_heatmap_super
        self.shred_modules_df_2 = df_subsetDEG

        if self.save_output:
            self.heatmap_matrix.to_csv(self.output_folder + "heatmap_matrix.txt", sep = "\t")
            self.heatmap_matrix_super.to_csv(self.output_folder + "heatmap_matrix_superbin.txt", sep = "\t")
            self.shred_modules_df_2.to_csv(self.output_folder + "genemodules_heatmap.txt", sep = "\t")

        # return [df_heatmap, df_heatmap_super, df_subsetDEG]


    def draw_heatmap(self):
        """
        Draw heatmap based on gene modules generated. Run after create_heatmap_matrix
        """
        try:
            self.heatmap_matrix
        except AttributeError:
            raise Exception("Heatmap should be created using create_heatmap_matrix prior to drawing it.")

        heatmap(self, bin_type = "bin", save_output = self.save_output)
        heatmap(self, bin_type = "superbin", save_output = self.save_output)


    def enrich_modules(self, 
                       categories = ["GeneOntologyMolecularFunction", 
                                    "GeneOntologyBiologicalProcess", 
                                    "GeneOntologyCellularComponent", 
                                    "Pathway", 
                                    "MousePheno"], 
                       ranked = False):
        """
        Do enrichment for gene modules generated. ToppGene curated enrichment knowledge is used.

        Parameters
        ----------
        categories
                Gene enrichment analysis categories used. 
                Current curation include GeneOntologyBiologicalProcess, GeneOntologyCellularComponent, 
                                         GeneOntologyMolecularFunction, MousePheno, Pathway. 
                More categories will be added later. For more information, please check https://toppgene.cchmc.org/
        ranked
                Whether use pre-ranked enrichment or not. Default is False. 
                For more information, please check package gseapy (https://gseapy.readthedocs.io/en/latest/introduction.html).
        """
        try:
            df_subsetDEG = self.shred_modules_df_2 
        except AttributeError:
            raise Exception("Gene modules should be generated using create_heatmap_matrix prior to gene module enrichment.")
        self.enrich_ranked = ranked
        
        df_enrich_output_all = pd.DataFrame()
        for status in np.unique(df_subsetDEG["Status"]):
            df_genelist = df_subsetDEG.loc[df_subsetDEG["Status"] == status, :]
            
            target = df_genelist.iloc[0, 0]
            reference = df_genelist.iloc[0, 7]
            module_name = (str(target) + "|" + str(reference)) if reference != "" else target

            df_genelist["genes"] = df_genelist.index.values # reformat the table
            df_genelist = df_genelist[["genes", "scores"]]
            
            # do enrichment
            if ranked == True:
                df_enrich_output = module_enrich_ranked(df_genelist, terms = categories)
            else:
                df_genelist = df_genelist[["genes"]]
                df_enrich_output = module_enrich(df_genelist, terms = categories)

            df_enrich_output["module"] = module_name 
            df_enrich_output_all = pd.concat([df_enrich_output_all, df_enrich_output], axis = 0)

        if self.save_output:
            self.df_module_enrichment = df_enrich_output_all
            self.df_module_enrichment.to_csv(self.output_folder + "enrichment.txt", sep = "\t")

        # return df_enrich_output_all

    
    def toppcluster(self, draw_plot = True):
        """
        run toppcluster for all modules and do clustering on heatmap. Check ToppCluster (https://toppcluster.cchmc.org/) for more details.

        Parameters
        ----------
        draw_plot
                Whether to draw a toppcluster map for comparative enrichment analysis.
        """
        try:
            self.df_module_enrichment
        except AttributeError:
            raise Exception("Module enrichment should be done prior to ToppCluster map generation.")
            
        df_toppcluster_modulemap = apply_toppcluster(self)
        df_toppcluster_modulemap.to_csv(self.output_folder + "enrichment_toppcluster.txt", sep = "\t")

        fig = sns.clustermap(df_toppcluster_modulemap, vmin = 0, vmax = 10, 
                                figsize=(15,10), xticklabels = False, yticklabels = True, cmap = "bwr")
        if draw_plot == True:
            fig.savefig(self.output_folder + "figures/toppcluster_map.png")

        # return df_toppcluster_modulemap


    def createJson(self):
        """
        create json file. It's currently under construction.
        """
        #self.json = json
        return 1
    

    def createGCT(self):
        """
        create GCT files based on both (super) heatmap matrix and bin / module metadata.
        """
        try:
            self.bin_metadata
            self.shred_modules_df_2
            self.heatmap_matrix
        except:
            raise Exception("Heatmap generated should be done prior to GCT file generation.")

        for bin_input in ["bin", "superbin"]:
            createGCT_file(self, bin_input)

    
    def toppcell_batchRun(self, 
                          top_n_genes = 200, 
                          heatmap_output_name = "heatmap.png",
                          enrich_categories = ["GeneOntologyMolecularFunction", 
                                                "GeneOntologyBiologicalProcess", 
                                                "GeneOntologyCellularComponent", 
                                                "Pathway", 
                                                "MousePheno"],
                          enrich_ranked = False,
                          toppcluster_run = True,
                          toppcluster_drawplot = True,
                          createGCT = True):
        """
        Run the whole pipeline in one code. Gene modules, heatmap and ToppCluster map can be generated in one run.

        Parameters
        ----------
        top_n_genes
                Number of genes selected for gene module generation.
        heatmap_output_name
                Output name of the heatmap.
        enrich_categories
                Categories used for enrichment.
        enrich_ranked
                Whether use pre-ranked enrichment or not. Default is False. 
        toppcluster_run
                Whether generate toppcluster map or not. Default is True.
        toppcluster_drawplot
                Whether draw toppcluster map or not. Default is True.
        createGCT
                Whether create GCT. Default is True
        """

        self.do_shredplan()
        self.create_heatmap_matrix(top_n_genes = top_n_genes)
        self.draw_heatmap(heatmap_output_name)
        self.enrich_modules(categories = enrich_categories, ranked = enrich_ranked)
        
        if toppcluster_run:
            self.toppcluster(draw_plot = toppcluster_draw_plot)
        if createGCT:
            self.createGCT()


def get_all_terms(terms):
    all_terms = []
    for term in terms:
        if "+" not in term:
            all_terms.append(term)
        else:
            all_terms += term.split("+")


def createGCT_file(shred, type_bin):
    """
    create GCT file based on both heatmap matrix and bin / module metadata
    """
    # initialization
    if type_bin == "bin":
        bin_meta = shred.bin_metadata
        heatmap_matrix = shred.heatmap_matrix
    elif type_bin == "supbin":
        bin_meta = shred.superbin_metadata
        heatmap_matrix = shred.heatmap_matrix_super
    else:
        raise Exception("type should be either bin or superbin")
    heatmap_matrix = shred.heatmap_matrix

    # format table
    bin_meta = bin_meta.loc[list(heatmap_matrix.columns), :]
    bin_meta["bin_id"] = bin_meta.index.values

    df_module.reset_index(level = 0, inplace = True)
    df_module.rename({"names": "Genes"}, axis = 1, inplace = True)

    df_new_heatmap = pd.DataFrame(data = heatmap_matrix.values,
                                    index = pd.MultiIndex.from_frame(df_module),
                                    columns = pd.MultiIndex.from_frame(bin_meta))

    gct_output_name = "heatmap_matrix_GCT.txt" if type_bin == "bin" else "heatmap_matrix_superbin_GCT.txt"
    df_new_heatmap.to_csv(shred.output_folder + gct_output_name, sep = "\t")
    
    # reformat the table into GCT3 type.
    second_line = str(heatmap_matrix.shape[0]) + "\t" + str(heat_matrix.shape[1]) + "\t" + str(df_module.shape[1]-1) + "\t" + str(bin_meta.shape[1]-1) + "\n"
    
    third_row = ""
    for i, col in enumerate(df_module.columns):
        if i != (df_module.shape[1] - 1):
            third_row += (col + "\t")
        else:
            third_row += (col + "\n")
    with open(shred.output_folder + gct_output_name, "r") as f:
        lines = f.readlines()
    with open(shred.output_folder + gct_output_name, "w") as f:
        f.write("#1.3\n")
        f.write(second_line)
        f.write(third_line)
        for line in lines:
            if not line.startswith("Genes"):
                f.write(line)