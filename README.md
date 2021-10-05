# ToppCell in Python
ToppCell (https://toppcell.cchmc.org/) is a web portal designed for biologists and bioinformaticians to explore single-cell RNA-seq data swiftly in the era of single cell and high throughput sequencing. See [Jin et al. 2021](https://www.sciencedirect.com/science/article/pii/S258900422101083X) for more details. This python version is still in development and your suggestions for improvement are really appreciated. 

### Main functions
- Creation of cell-label-based hierarchical gene modules with most significant signatures for single-cell data.
- Visualization of gene modules in comprehensive way.
- Large scale gene enrichment analysis supported by ToppGene.
- ToppCluster-based comparable gene enrichment analysis.

### Reference
```
Jin, Kang, et al. "An Interactive Single Cell Web Portal Identifies Gene and Cell Networks in COVID-19 Host Responses." Iscience (2021): 103115.
```

### Installation
```python
git clone https://github.com/KANG-BIOINFO/ToppCellPy.git
cd ToppCellPy
pip install .
```

### Tutorial
- [Example on COVID-19 data]()
- [Example on IFN-stimulated PBMC data](https://nbviewer.jupyter.org/github/KANG-BIOINFO/ToppCellPy/blob/main/test/IFN-stimulated%20PBMC%20example.ipynb)

### Quick Start
```python
import scanpy as sc
import numpy as np
import ToppCellPy as tp
```

Load example data
```python
adata = sc.read("batch2_all_normalized_filtered.h5ad")
```

Set up shred object run shred plan
```python
shred = tp.Shred(adata = adata,
            shred_plan = ["stim", "cell", "stim+cell|stim"],
            bin_group = ["stim", "cell"],
            order_bins = None,
            order_modules = None,
            method = "wilcoxon")
shred.do_shredplan()
```

Create heatmap data frame and draw heatmap
```python
shred.create_heatmap_matrix()
shred.draw_heatmap()
```

Do enrichment for all modules
```python
shred.enrich_modules(categories = ["GeneOntologyCellularComponent"])
```

Draw ToppCluster plot
```python
shred.toppcluster()
```

An alternative option for running the whole pipeline above is using the one-stop code below:
```python
shred = tp.Shred(adata = adata,
            shred_plan = ["stim", "cell", "stim+cell|stim"],
            bin_group = ["stim", "cell"],
            order_bins = None,
            order_modules = None,
            method = "wilcoxon")
shred.toppcell_batchRun(top_n_genes = 200,
                        enrich_categories = ["GeneOntologyCellularComponent"],
                        enrich_ranked = False,
                        toppcluster_run = True,
                        createGCT = True)
```