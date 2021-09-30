# ToppCell in Python
ToppCell (https://toppcell.cchmc.org/) is a web portal designed for biologists and bioinformaticians to explore single-cell RNA-seq data swiftly in the era of single cell and high throughput sequencing. Now we enable the visualization function using Python, making it flexible and standalone. It's still improving and your suggestions are really appreciated.

### Main function
- Creation of cell-label-based hierarchical gene modules with most signatures for single-cell data.
- Automated cell type annotation with transfer learning.
- Visualization of gene modules in comprehensive way.
- Large scale gene enrichment analysis supported by ToppGene.
- ToppCluster-based comparable gene enrichment analysis.
- Others...

### Installation
```python
git clone https://github.com/KANG-BIOINFO/ToppCell_Python.git
cd ToppCell_Python
pip install .
```

### Tutorial
- Example on COVID-19 data

### Quick Start
```python
import scanpy as sc
import numpy as np
import ToppCell_Python as tp
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
shred.draw_heatmap("heatmap.png")
```

Do enrichment for all modules
```python
shred.enrich_modules(categories = ["GeneOntologyCellularComponent"])
```

Draw ToppCluster plot
```python
shred.toppcluster()
```


