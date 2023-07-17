import os
from scipy.io import mmwrite
from scipy.spatial.distance import pdist
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import cellrank as cr
import scvelo as scv
import igraph as ig

import cmocean

os.chdir('/home/jfleck/projects/early/')

#### Get velocity ####
rna = sc.read('data/RNA_ATAC/subsets/RNA_all.h5ad')
velo_all = sc.read('/local1/USERS/jfleck/data/EARLY/RNA_ATAC/all_rnatac_velo.h5ad')
genes_use = rna.var.index.intersection(velo_all.var.index)
velo_all = velo_all[cells_use, genes_use]
rna = rna[:, genes_use]

rna.var = velo_all.var
rna.layers['spliced'] = velo_all.layers['spliced']
rna.layers['unspliced'] = velo_all.layers['unspliced']

rna_var_df = sc.pp.highly_variable_genes(rna, n_top_genes=2000, inplace=False)
rna_highvar_genes = rna.var.index[rna_var_df.highly_variable]

var_genes_use = rna_highvar_genes.difference(cc_genes)

rna_var = rna[:, var_genes_use]

#### Neighbors and moments based on CSS integration ####
ncss = rna.obsm['X_css'].shape[1]
sc.pp.neighbors(rna_var, n_neighbors=30, n_pcs=ncss, use_rep='X_css')
scv.pp.moments(rna_var, n_pcs=ncss, n_neighbors=30, use_rep='X_css')
scv.tl.velocity(rna_var, mode='stochastic')
scv.tl.velocity_graph(rna_var)

rna_var.write(f'data/RNA_ATAC/subsets/all_velo_stoch_2000f.h5ad')


#### Cellrank prep ####
velo_full = sc.read('data/RNA_ATAC/velocity/all_velo_stoch_2000f.h5ad')
# velo_full.write('data/RNA_ATAC/subsets/all_velo_stoch_2000f.h5ad')

velo_full.obs['neuron_type'] = rna.obs['neuron_type']
velo_full.obs['RNA_snn_res.10'] = rna.obs['RNA_snn_res.10']
velo_full.obs['RNA_snn_res.20'] = rna.obs['RNA_snn_res.20']
velo_full.obs['RNA_snn_res.2'] = rna.obs['RNA_snn_res.2']
velo_full.obs['velocity_pseudotime'] = rna.obs['velocity_pseudotime']
velo_full.obs['pseudotime_ranks'] = rna.obs['pseudotime_ranks'] / rna.obs['pseudotime_ranks'].max()

sc.pl.scatter(velo_full, basis='umap', color=['RNA_snn_res.2', 'neuron_type'], legend_loc='on data')

age_colors = rev(c(
    '#264653','#287271','#2a9d8f','#5aa786','#8ab17d',
    '#A2B679', '#babb74','#e9c46a','#efb366','#ee8959','#e76f51'))

velo_full.obs['age'] = velo_full.obs['age'].astype(int).astype(str).astype('category')

scv.pl.velocity_embedding_grid(velo_full, basis='umap', color='RNA_snn_res.2', legend_loc='on data')


#### Get velocity pseudotime ####
root_cluster = velo_full[velo_full.obs['RNA_snn_res.2']=='23']
root_cell = root_cluster[root_cluster.obsm['X_umap'][:,1].argmax(), :].obs.index.values[0]
root_idx = np.where(velo_full.obs.index==root_cell)[0][0]
velo_full.obs['is_root'] = velo_full.obs.index==root_cell

scv.tl.velocity_pseudotime(velo_full, root_key=root_idx)

sc.pl.scatter(velo_full, basis='umap', color=['velocity_pseudotime'], legend_loc='on data')

velocity_meta = velo_full.obs[['velocity_self_transition', 'velocity_pseudotime']]
velocity_meta['pseudotime_ranks'] = velocity_meta['velocity_pseudotime'].rank().astype(int)
velocity_meta.to_csv('data/RNA_ATAC/RNA_ATAC_full_velocity_pt.tsv', sep='\t')

# velo_full.write('data/RNA_ATAC/velocity/all_velo_stoch_2000f.h5ad')
# velo_full = sc.read('data/RNA_ATAC/velocity/all_velo_stoch_2000f.h5ad')

#### Get transition probabilities to neuronal states with cellrank ####
# sc.pp.neighbors(velo_full, n_neighbors=30, use_rep='X_css')

vk = cr.tl.kernels.VelocityKernel(velo_full)
vk.compute_transition_matrix()

pk = cr.tl.kernels.PseudotimeKernel(velo_full, time_key = 'velocity_pseudotime')
pk.compute_transition_matrix(threshold_scheme='hard')

ck = cr.tl.kernels.ConnectivityKernel(velo_full)
ck.compute_transition_matrix()

combined_kernel = 0.5 * vk + 0.5 * ck

g = cr.tl.estimators.GPCCA(combined_kernel, write_to_adata=False)

terminal_states = velo_full.obs['neuron_type'].astype(str)
terminal_states[terminal_states=='NA'] = np.nan
terminal_states = terminal_states.astype('category')

g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities(n_jobs=16)
g.plot_absorption_probabilities(mode='embedding', same_plot=False, basis='umap')

velo_full.obs['to_ctx'] = velo_full.obsm['to_terminal_states']['ctx_ex'].X
velo_full.obs['to_nt'] = velo_full.obsm['to_terminal_states']['mesen_ex'].X
velo_full.obs['to_ge'] = velo_full.obsm['to_terminal_states']['ge_in'].X

velo_full.obs['to_ctx_ranks'] = velo_full.obs['to_ctx'].rank() 
velo_full.obs['to_nt_ranks'] = velo_full.obs['to_nt'].rank() 
velo_full.obs['to_ge_ranks'] = velo_full.obs['to_ge'].rank() 

velo_full.obs['to_ctx_ranks'] = velo_full.obs['to_ctx_ranks'] / max(velo_full.obs['to_ctx_ranks'])
velo_full.obs['to_nt_ranks'] = velo_full.obs['to_nt_ranks'] / max(velo_full.obs['to_nt_ranks'])
velo_full.obs['to_ge_ranks'] = velo_full.obs['to_ge_ranks'] / max(velo_full.obs['to_ge_ranks'])

sc.pl.scatter(velo_full, basis='umap', color='to_nt_ranks')
sc.pl.scatter(velo_full, basis='umap', color=['to_ge', 'to_ctx', 'to_nt'])
sc.pl.scatter(velo_full, basis='umap', color=['to_ge_ranks', 'to_ctx_ranks', 'to_nt_ranks'])

# velo_full.write('data/RNA_ATAC/subsets/_tmp_all_cr_8vk_2pk_0ck.h5ad')
# velo_full = sc.read('data/RNA_ATAC/subsets/_tmp_all_cr_8vk_2pk_0ck.h5ad')

cellrank_meta = velo_full.obs[['to_ctx', 'to_nt', 'to_ge', 'to_ge_ranks', 'to_ctx_ranks', 'to_nt_ranks', 'pseudotime_ranks', 'terminal_states_probs']]
cellrank_meta.to_csv('data/RNA_ATAC/velocity/RNA_ATAC_full_cellrank_probs.tsv', sep='\t')

trans_mat = combined_kernel._transition_matrix
mmwrite('data/RNA_ATAC/velocity/velo_cr_transition.mtx', trans_mat)


#### PAGA with trans probs ####
root_cluster = (velo_full.obs['RNA_snn_res.20']=='59').astype(float)
velo_full.obs['root_prob'] = root_cluster
velo_full.obs['RNA_snn_res.20'] = velo_full.obs['RNA_snn_res.20'].astype('category')

sc.pl.scatter(velo_full, basis='umap', color=['root_prob'])

scv.tl.paga(
    velo_full,
    groups='RNA_snn_res.20',
    use_time_prior='velocity_pseudotime'
)


scv.pl.paga(velo_full, basis='umap')
con_mat = velo_full.uns['paga']['connectivities']

mmwrite('data/RNA_ATAC/velocity/velo_paga_connectivities.mtx', con_mat)















