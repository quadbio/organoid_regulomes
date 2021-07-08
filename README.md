## Code repo for the manuscript: "Inferring and perturbing cell fate regulomes in human cerebral organoids"

The repo is structured as follows:

* `integration/` contains scripts and functions used to integrate the RNA-seq and ATAC-seq data
	- `get_bipartite_matches.py` finds matching cells given a h5ad file with a common embedding (in our case CCA) between the datasets.
	- `matching.py` contains helper functions for the matching script
* `pando/` contains scripts and functions used to infer the gene regulatory network
	- `atac.R` contains functions to manipulate and plot ATAC-seq data
	- `pseudocells.R` contains the functions to construct pseudocells as used in [Kanton et al. 2019](https://www.nature.com/articles/s41586-019-1654-9)
	- `fit_multiome_glm.R` contains a script to do model fitting for GRN inference. This is now also implemented in the R package [Pando](https://github.com/quadbiolab/Pando)
	- `models.R` contains helper functions for model fitting.
* `crop_seq/` contains scripts and functions used for the analysis of the CROP-seq data
	- `get_guide_umis.py` is a script to extract guide UMIs with high read support given a cellranger output (`molecule_info.h5`) using a GMM-based approach
	- `ko_inference.R` performs estimation of perturbation probabilities using a linear model-based approach as described in the [Perturb-seq paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5181115/). It's essentially a R implementation of the approach used in [MIMOSCA](https://github.com/asncd/MIMOSCA).
	- `enrichment.R` implements statistical tests for detection of composition changes
	- `plots.R` has some functions to plot CROP/Perturb-seq data.
* `utils/` contaions general utility functions
	- `de.R` contains functions to perform differential expression and differential accessibility analysis
	- `utils.R` has some other useful functions
