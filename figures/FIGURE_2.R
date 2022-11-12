source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

setwd('~/projects/early/')

library(ggraph)
library(tidygraph)
library(tidytree)
library(fastmatch)
library(biomaRt)
library(ggh4x)
library(doParallel)
library(IRanges)
library(GenomicRanges)

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarize

registerDoParallel(10) 

#### Load data ####
rnatac <- read_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2.1_srt.rds')
# rnatac %>% write_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2.1_srt.rds')

cluster_summaries <- read_rds('data/RNA_ATAC/RNA_ATAC_res20_summaries.rds')
feature_sets <- read_rds('data/gene_sets/RNA_feature_sets.rds')


gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')
gene_annot <- gene_annot[gene_annot$gene_name %in% feature_sets$grn_features]
peak_ranges <- StringToGRanges(rownames(rnatac@assays$peaks))

atac <- read_rds('data/ATAC/ATAC_all_merged/ATAC_all_merged_v3.1_srt.rds')
# atac %>% write_rds('data/ATAC/ATAC_all_merged/ATAC_all_merged_v3.1_srt.rds')


## Peak matches
peak_subsets <- read_rds('data/RNA_ATAC/RNA_ATAC_peak_subsets_v2.rds')
grn_regions <- peak_subsets$cons_cre_noex
grnreg2peaks <- rownames(rnatac@assays$peaks@data)[grn_regions$peaks]
names(grnreg2peaks) <- GRangesToString(grn_regions$ranges)


## Select allowed genes 
mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl')
all_coding_genes <- getBM(attributes = c( 'hgnc_symbol'), filters = c('biotype'), values = list(biotype='protein_coding'), mart = mart)
good_genes <- all_coding_genes$hgnc_symbol
good_genes <- good_genes[!str_detect(good_genes, '^MT-|^RP|^HIST')]

motif2tf <- read_tsv('data/tf_matches/motif2tf_early_all.tsv') 
# motif2tf %>% write_tsv('data/tf_matches/motif2tf_early_all.tsv') 

lineage_graph <- read_rds('data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')
hc_meta <- read_tsv('data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv') %>% 
    mutate(cluster=as.character(cluster))
cluster_meta <- as_tibble(lineage_graph)

umap_meta <- Reductions(rnatac, slot='umap')@cell.embeddings %>% 
    as_tibble(rownames='cell')

meta <- rnatac@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(umap_meta)

## Gene and GRN stuff
modules <- read_tsv('data/RNA_ATAC/grn/networks/v5_d100k_cons_cre_noex_modules.tsv')
gene_cor <- read_tsv('data/RNA_ATAC/grn/networks/coexpression_corr.tsv')
gene_meta <- read_tsv('data/gene_sets/gene_scores.tsv')
tf_graph <- read_rds('data/RNA_ATAC/grn/tf_graph.rds')
tf_graph_meta <- tf_graph %N>% as_tibble()



plot_df <- tf_graph_meta %>% 
    mutate(pt_bins=cut(mean_pt, c(0,0.35,0.5,0.8,1), labels=c('psc/nect', 'nepi', 'npc', 'neuron')))

p1 <- ggplot(plot_df, aes(pt_bins, central_pr)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(x='stage', y='pagerank centrality')

p2 <- ggplot(plot_df, aes(pt_bins, central_eig)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(x='stage', y='eigen centrality')

p3 <- ggplot(plot_df, aes(pt_bins, ngenes)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(x='stage', y='module size (outdegree)')

p1 | p2 | p3
ggsave('plots/RNA_ATAC/grn/centrality_stages_boxplot.png')



## DE results
module_score_second_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_pos_second_split.tsv')
module_pos_peaks_second_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_pos_peaks_second_split.tsv')
module_neg_second_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_neg_second_split.tsv')

module_score_first_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_pos_first_split.tsv')
module_pos_peaks_first_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_pos_peaks_first_split.tsv')
module_neg_first_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_neg_first_split.tsv')

module_score_lin_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_pos_lineage.tsv')
module_pos_peaks_lin_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_pos_peaks_lineage.tsv')
module_neg_lin_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_module_neg_lineage.tsv')



#### Get gene and peak modules ####
module_pos <- modules %>% 
    filter(estimate>0) %>% 
    group_by(tf) %>% filter(n()>5) %>% 
    group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>% 
    map(function(x) x$target)

module_neg <- modules %>% 
    filter(estimate<0) %>% 
    group_by(tf) %>% filter(n()>5) %>% 
    group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>% 
    map(function(x) x$target)

regions_pos <- modules %>% 
    filter(estimate>0) %>% 
    group_by(tf) %>% filter(n()>5) %>% 
    group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>% 
    map(function(x) unlist(str_split(x$peaks, ';')) %>% str_replace_all('_', '-'))

regions_neg <- modules %>% 
    filter(estimate<0) %>% 
    group_by(tf) %>% filter(n()>5) %>% 
    group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>% 
    map(function(x) unlist(str_split(x$peaks, ';')) %>% str_replace_all('_', '-'))

peaks_pos <- regions_pos %>% map(function(x) unique(grnreg2peaks[x]))
peaks_neg <- regions_neg %>% map(function(x) unique(grnreg2peaks[x]))

all_modules <- list(
  'genes_pos' = module_pos,
  'genes_neg' = module_neg,
  'regions_pos' = regions_pos,
  'regions_neg' = regions_neg,
  'peaks_pos' = peaks_pos,
  'peaks_neg' = peaks_neg
)

all_modules %>% write_rds('data/RNA_ATAC/grn/networks/tf_peak_gene_modules.rds')


#### Network viz ####
ggraph(tf_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=mean_ranked_pt, size=central_pr), shape=21, stroke=0.2) +
    # geom_node_text(aes(label=name, filter=name%in%genes_show, size=central_pr), repel=T) +
    # geom_node_text(aes(label=name), size=2) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void() + no_legend()
ggsave('plots/paper/fig2/tf_network_pt_pr_umap.png', width=15, height=15, units='cm')



ggraph(tf_graph, x=UMAP1, y=UMAP2) + 
  geom_edge_diagonal(color='gray', width=0.1, alpha=0.1) + 
  geom_node_point(color='darkgrey', size=2) +
  theme_void()
ggsave('plots/paper/fig2/tf_network_minimal_umap.png', width=15, height=15, units='cm')


ggraph(tf_graph, x=UMAP1, y=UMAP2) + 
  # geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) + 
  scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
  scale_edge_alpha_continuous(range=c(0.01,0.4), limits=c(2,20)) +
  geom_node_point(aes(fill=mean_ranked_pt, size=central_pr), shape=21, stroke=0.2) +
  geom_node_text(aes(label=name), size=2, repel=T, max.overlaps=1000) +
  scale_fill_viridis(option='magma', direction = -1) +
  theme_void() + no_legend()
ggsave('plots/paper/fig2/tf_network_all_labels_umap.pdf', width=25, height=25, units='cm')


ggraph(tf_graph, x=UMAP1, y=UMAP2) + 
    # geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) + 
    geom_node_point(size=1, alpha=1) +
    geom_node_text(aes(label=name, filter=name%in%genes_show), size=5/ggplot2::.pt, repel=T) +
    theme_void()
ggsave('plots/paper/fig2/tf_network_labels.pdf', width=7.8, height=7.8, units='cm')



ggraph(arrange(tf_graph, ngenes), x=UMAP1, y=UMAP2) + 
    geom_node_point(aes(fill=ngenes), shape=21, stroke=0.2, size=2.5) +
    # geom_node_text(aes(label=name, filter=name%in%genes_show), size=2, repel=T) +
    scale_fill_gradientn(colors=rev(pals::ocean.gray(100)), limits=c(1,800)) +
    scale_size_continuous(range=c(1,5)) +
    theme_void() 
ggsave('plots/paper/fig2/tf_network_ngenes_umap.png', width=8, height=8, units='cm')


#### Different centrality metrics ####
module_graph <- as_tbl_graph(modules) %N>% 
  mutate(
    # central_katz  = centrality_katz(),
    central_close  = centrality_closeness(),
    central_alpha = centrality_alpha(),
    central_exp = centrality_expected(),
    central_inf = centrality_information(),
    central_decay = centrality_decay()
  )

modules <- modules %>% 
  select(-corr) %>% 
  left_join(gene_cor, by=c('tf'='source', 'target'))


tf_meta <- modules %>% ungroup() %>% 
  distinct('name'=tf, ngenes)
target_meta <- modules %>% ungroup() %>% 
  distinct('name'=target, ntfs)

grn_graph <- as_tbl_graph(modules) %>% 
  activate(edges) %>% 
  mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>% 
  activate(nodes) %>% left_join(tf_meta) %>% left_join(target_meta) %>% 
  inner_join(select(as_tibble(activate(tf_graph, nodes)), name, UMAP1, UMAP2)) %>% 
  mutate(
    central_pr = centrality_pagerank(weights = estimate),
    central_betw = centrality_betweenness(),
    central_eig = centrality_eigen(),
    centrality_auth = centrality_authority(),
    central_hub = centrality_hub(weights=estimate),
    central_outdegree = centrality_degree(mode='out'),
    central_indegree = centrality_degree(mode='in'),
    central_harm = centrality_closeness_harmonic(),
    central_gen = centrality_closeness_generalised(alpha=0.1),
    central_int = centrality_integration(),
    central_comm = as.double(centrality_communicability()),
    central_curr = centrality_betweenness_current()
  ) %>% to_undirected() %>% 
  mutate(
    group_louvain=group_louvain()
  ) 

centr_metrics <- colnames(as_tibble(grn_graph)) %>% {.[str_detect(., 'central')]}
plots <- map(centr_metrics, function(metric){
  ggraph(grn_graph, x=UMAP1, y=UMAP2) + 
    geom_node_point(aes_string(fill=metric, size=metric), shape=21, stroke=0.2) +
    scale_fill_viridis(option='viridis') +
    theme_void() + no_legend() + ggtitle(metric)
})
wrap_plots(plots)
ggsave('plots/paper/fig2/tf_network_centrality_umap.png', width=20, height=20)



#### Annotate modules based on activity ####
rnatac@active.assay <- 'RNA'
rnatac$this <- rnatac$pseudotime_ranks > 0.4
feature_plot(rnatac, features=c('S.Score', 'G2M.Score', 'VIM', 'is_npc'))
dim_plot(rnatac, group.by=c('is_mesoderm'))

rnatac$is_neuron <- !is.na(rnatac$neuron_type)
rnatac$is_npc <- (rnatac$pseudotime_ranks > 0.4) & (is.na(rnatac$neuron_type))
rnatac$diff_state <- case_when(
  rnatac$is_neuron ~ 'neuron',
  rnatac$is_npc ~ 'npc',
  rnatac$is_mesoderm ~ 'mesenchyme',
  TRUE ~ 'early'
)


### Mean module score
cc_module_summary <- aggregate_matrix(t(rnatac[['module_score']]@data), groups=rnatac$Phase)
neuron_type_module_summary <- aggregate_matrix(t(rnatac[['module_score']]@data), groups=rnatac$neuron_type)
neuron_module_summary <- aggregate_matrix(t(rnatac[['module_score']]@data), groups=rnatac$diff_state)
age_module_summary <- aggregate_matrix(t(rnatac[['module_score']]@data), groups=rnatac$age) %>% 
  {rownames(.) <- paste0('age_', rownames(.)); .}

# cc_module_summary <- aggregate_matrix(t(rnatac[['RNA']]@data[rownames(rnatac[['module_score']]), ]), groups=rnatac$Phase)
# neuron_type_module_summary <- aggregate_matrix(t(rnatac[['RNA']]@data[rownames(rnatac[['module_score']]), ]), groups=rnatac$neuron_type)
# neuron_module_summary <- aggregate_matrix(t(rnatac[['RNA']]@data[rownames(rnatac[['module_score']]), ]), groups=rnatac$diff_state)

cc_module_df <- cc_module_summary %>% t() %>% 
  as_tibble(rownames='name') 

ntype_module_df <- neuron_type_module_summary %>% t() %>% {.[, colnames(.)!='NA']} %>% 
  as_tibble(rownames='name') 

diff_module_df <- neuron_module_summary %>% t() %>% 
  as_tibble(rownames='name')

age_module_df <- age_module_summary %>% t() %>% 
  as_tibble(rownames='name') 

module_meta <- purrr::reduce(list(cc_module_df, ntype_module_df, diff_module_df, age_module_df), inner_join)
  
grn_df <- grn_graph %N>% 
  inner_join(module_meta) %>% 
  as_tibble()

summ_vars <- colnames(module_meta)[-1]
plots <- map(summ_vars, function(metric){
  clim <- max(abs(grn_df[[metric]]))
  ggplot(arrange_at(grn_df, metric), aes(x=UMAP1, y=UMAP2)) + 
    geom_point(aes_string(fill=metric, size=metric), shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=grad(pals::brewer.greys, 0.5), limits=c(-clim,clim)) +
    scale_size_continuous(range=c(0.1,3)) +
    theme_void() + no_legend() + ggtitle(metric)
})
wrap_plots(plots)

ggsave('plots/paper/fig2/tf_network_centrality_umap.png', width=20, height=20)


### DE FC 
cc_module_de <- de(rnatac, assay='module_score', groups='Phase')
neuron_module_de <- de(rnatac, assay='module_score', groups='neuron_type')
diff_module_de <- de(rnatac, assay='module_score', groups='diff_state')

module_de_meta <- map(list(cc_module_de, neuron_module_de, diff_module_de), select, 'name'=feature, group, fc) %>% 
  bind_rows() %>% pivot_wider(names_from = 'group', values_from='fc')

grn_df <- grn_graph %N>% 
  inner_join(module_meta) %>% 
  as_tibble()

summ_vars <- colnames(module_de_meta)[-1]
plots <- map(summ_vars, function(metric){
  clim <- max(abs(grn_df[[metric]]))
  ggplot(arrange_at(grn_df, metric), aes(x=UMAP1, y=UMAP2)) + 
    geom_point(aes_string(fill=metric, size=metric), shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=grad(pals::brewer.greys, 0.5)) +
    scale_size_continuous(range=c(0.1,3)) +
    theme_void() + no_legend() + ggtitle(metric)
})
wrap_plots(plots)



modules %>% filter(target=='GLI3')




#### Graph branches ####
plot_graph <- lineage_graph %>% 
    mutate(lineage_show=case_when(
        lineage%in%c('ctx', 'ge') ~ 'telencephalon', 
        lineage=='early' ~ 'other',
        T ~ lineage
    ))

ggraph(plot_graph, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.2, width=0.5) +
    geom_node_point(aes(fill=lineage_show), size=3, shape=21, stroke=0.2) +
    scale_fill_manual(values=lineage_colors) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend()
ggsave('plots/paper/fig2/graph_first_split.png', width=4.3, height=3.7)


plot_graph <- lineage_graph %>% 
    mutate(lineage_show=ifelse(lineage%in%c('ctx', 'ge') & pseudotime_ranks>0.6, lineage, 'other'))

ggraph(plot_graph, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.2, width=0.5) +
    geom_node_point(aes(fill=lineage_show), size=3, shape=21, stroke=0.2) +
    scale_fill_manual(values=lineage_colors) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend()
ggsave('plots/paper/fig2/graph_second_split.png', width=4.3, height=3.7)


ggraph(plot_graph, x=FR1, y=FR2) +
  geom_edge_link(alpha=0.2, width=0.5) +
  geom_node_point(aes(fill=pseudotime_ranks>0.3), size=3, shape=21, stroke=0.2) +
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_void() + no_legend()



#### Gene vs Peak FC ####
#### First split ####
score_first <- module_score_first_de %>% 
    filter(group%in%c('telencephalon')) 

peaks_first <- module_pos_peaks_first_de %>% 
    filter(group%in%c('telencephalon'))

label_genes <- c('OTX1', 'SP1', 'LMX1A', 'SP5', 'LMX1B', 'TEAD1', 'PAX8', 'FEZF1', 'GSX2', 'DLX2', 'DLX1', 'SATB1', 'MEIS3', 'ZFHX3', 'FOXG1')

module_meta <- modules %>% 
    group_by(tf) %>% 
    summarize(
        ngenes=ngenes[1],
        npeaks=length(unlist(str_split(peaks, ';')))
    ) %>% rename('gene'='tf')

score_vs_peaks <- inner_join(
    select(score_first, feature, fc, padj, auc),
    select(peaks_first, feature, fc, padj, auc),
    suffix = c('_genes', '_peaks'), by=c('feature')
) %>% inner_join(module_meta, by=c('feature'='gene')) %>% 
    mutate(label=feature%in%label_genes | abs(fc_genes)>0.05 & abs(fc_peaks)>0.2, both_sig=padj_peaks<1e-5 & padj_genes<1e-5) %>% arrange(!both_sig)

table(sign(filter(score_vs_peaks, both_sig)$fc_genes), sign(filter(score_vs_peaks, both_sig)$fc_peaks))
summary(lm(filter(score_vs_peaks, both_sig)$fc_genes~filter(score_vs_peaks, both_sig)$fc_peaks))

ggplot(score_vs_peaks, aes(fc_genes, fc_peaks, label=feature, fill=both_sig, alpha=both_sig)) +
    geom_hline(yintercept=0, color='darkgray', size=0.2) +
    geom_vline(xintercept=0, color='darkgray', size=0.2) +
    scale_fill_manual(values=c('grey', '#34495e')) +
    geom_point(shape=21, size=0.5, stroke=0.1, color='white') +
    geom_rangeframe(data=tibble(x=c(-0.3, 0.3), y=c(-0.6, 0.6)), mapping=aes(x=x, y=y), inherit.aes=F) +
    scale_x_continuous(limits = c(-0.35, 0.35), breaks=seq(-0.3, 0.3, 0.3)) +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks=seq(-0.6, 0.6, 0.3)) +
    scale_alpha_manual(values=c(0.5,1)) +
    # geom_text_repel(data=filter(score_vs_peaks, label), size=4) +
    no_legend() +
    article_text() +
    theme(
        panel.border = element_blank(),
        axis.ticks = element_line(size = 0.5)
    ) +
    labs(x='Module score diff.', y='Peak enrichment diff.')

ggsave('plots/paper/fig2/module_peaks_vs_genes_de_first_scatter.pdf', width=4, height=3.3, unit='cm')


ggplot(score_vs_peaks, aes(fc_genes, fc_peaks, label=feature, fill=both_sig, alpha=both_sig)) +
    geom_hline(yintercept=0, color='darkgray', size=0.2) +
    geom_vline(xintercept=0, color='darkgray', size=0.2) +
    scale_fill_manual(values=c('#d6dbdf', '#34495e')) +
    scale_x_continuous(limits = c(-0.3, 0.3), breaks=seq(-0.3, 0.3, 0.3)) +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks=seq(-0.6, 0.6, 0.3)) +
    scale_alpha_manual(values=c(0.5,1)) +
    geom_text(data=filter(score_vs_peaks, label), size=0.5) +
    theme_void() +
    no_legend() 


ggsave('plots/paper/fig2/module_peaks_vs_genes_de_first_scatter_labels.pdf', width=4, height=3.3, unit='cm')



#### Second split ####
rnatac@active.assay <- 'module_score'
rnatac@active.assay <- 'module_peaks_chromvar'
dim_plot(rnatac, group.by='lineage')
feature_plot(rnatac, features=c('ISL1', 'DLX5', 'SP8', 'L3MBTL3'), order=T)

rnatac$this <- rnatac$pseudotime_ranks>0.7
dim_plot(rnatac, group.by='this')


score_second <- module_score_second_de %>% 
    filter(group%in%c('ctx')) 

peaks_second <- module_pos_peaks_second_de %>% 
    filter(group%in%c('ctx'))

label_genes <- c('NFIB', 'EMX1', 'NR1D1', 'NEUROD2', 'DLX5', 'DLX2', 'ISL1', 'DMRT3', 'FOXN4', 'HEY1', 'GLI3', 'NKX2-1')

score_vs_peaks <- inner_join(
    select(score_second, feature, fc, padj, auc),
    select(peaks_second, feature, fc, padj, auc),
    suffix = c('_genes', '_peaks'), by=c('feature')
) %>% inner_join(module_meta, by=c('feature'='gene')) %>% 
    mutate(label=feature%in%label_genes | abs(fc_genes)>0.05 & (fc_peaks>0.8 | fc_peaks<(-0.4)), both_sig=padj_peaks<1e-5 & padj_genes<1e-5) %>% arrange(!both_sig)

table(sign(filter(score_vs_peaks, both_sig)$fc_genes), sign(filter(score_vs_peaks, both_sig)$fc_peaks))
summary(lm(filter(score_vs_peaks, both_sig)$fc_genes~filter(score_vs_peaks, both_sig)$fc_peaks))

ggplot(score_vs_peaks, aes(fc_genes, fc_peaks, label=feature, fill=both_sig, alpha=both_sig)) +
    geom_hline(yintercept=0, color='darkgray', size=0.2) +
    geom_vline(xintercept=0, color='darkgray', size=0.2) +
    scale_fill_manual(values=c('grey', '#34495e')) +
    geom_point(shape=21, size=0.5, stroke=0.1, color='white') +
    geom_rangeframe(data=tibble(x=c(-0.6, 0.3), y=c(-1, 2)), mapping=aes(x=x, y=y), inherit.aes=F) +
    scale_x_continuous(limits = c(-0.3, 0.3), breaks=seq(-0.3, 0.3, 0.3), na.value=-0.3) +
    scale_y_continuous(limits = c(-1.2, 1.5), breaks=seq(-1, 1.5, 0.5), na.value=1.5) +
    scale_alpha_manual(values=c(0.5,1)) +
    # geom_text_repel(data=filter(score_vs_peaks, label), size=4, max.overlaps=1000) +
    no_legend() +
    article_text() +
    theme(
        panel.border = element_blank(),
        axis.ticks = element_line(size = 0.5)
    ) +
    labs(x='Module score diff.', y='Peak enrichment diff.')

ggsave('plots/paper/fig2/module_peaks_vs_genes_de_second_scatter.pdf', width=3.6, height=3.3, unit='cm')


ggplot(score_vs_peaks, aes(fc_genes, fc_peaks, label=feature, fill=both_sig, alpha=both_sig)) +
    geom_hline(yintercept=0, color='darkgray', size=0.2) +
    geom_vline(xintercept=0, color='darkgray', size=0.2) +
    scale_fill_manual(values=c('#d6dbdf', '#34495e')) +
    # geom_point(shape=21, size=0.5, stroke=0.1) +
    scale_x_continuous(limits = c(-0.3, 0.3), breaks=seq(-0.3, 0.3, 0.3), na.value=-0.3) +
    scale_y_continuous(limits = c(-1.2, 1.5), breaks=seq(-1, 1.5, 0.5), na.value=1.5) +
    scale_alpha_manual(values=c(0.5,1)) +
    geom_text(data=filter(score_vs_peaks, label), size=0.5) +
    theme_void() +
    no_legend() 


ggsave('plots/paper/fig2/module_peaks_vs_genes_de_second_scatter_labels.pdf', width=3.6, height=3.3, unit='cm')



#### Fit LM to get R2 ####
rnatac$this <- rnatac$pseudotime_ranks > 0.4
dim_plot(rnatac, group.by=c('this', 'lineage', 'stage_manual_hclust'))

rnatac@active.assay <- 'module_score'
feature_plot(rnatac, features=c('SP9', 'HIF3A', 'DLX5', 'NR3C1', 'ADNP2'), order=T)


cells_use <- rnatac$pseudotime_ranks > 0.4 & rnatac$lineage%in%c('ctx', 'ge', 'nt')
meta <- rnatac@meta.data[cells_use, ] %>% 
  as_tibble(rownames='cell') %>% 
  group_by(lineage) %>% 
  mutate(pseudotime_align=rank(pseudotime_ranks)/n())
rnatac$pseudotime_align <- NA
rnatac$pseudotime_align[cells_use] <- meta$pseudotime_align



model_frame <- data.frame(
  lineage=factor(rnatac$lineage[cells_use]),
  pt=rnatac$pseudotime_align[cells_use]
)

module_mat <- t(rnatac[['module_score']]@data[, cells_use])

# module_mat <- t(rnatac[['RNA']]@data[rownames(rnatac[['module_score']]), cells_use])
colnames(module_mat) <- str_replace(colnames(module_mat), '-', '_')

module_fit_list <- pblapply(1:ncol(module_mat), function(i){
  model_formula <- reformulate(colnames(model_frame), response=colnames(module_mat)[i])
  model_data <- cbind(module_mat[,i,drop=F], model_frame)
  return(fit_glm(formula=model_formula, data=model_data))
})
names(module_fit_list) <- colnames(module_mat)
comb_model_coefs <- map_dfr(module_fit_list, function(x) x$coefs, .id='module')
comb_model_gof <- map_dfr(module_fit_list, function(x) x$gof, .id='module')

module_lin_fit_list <- pblapply(1:ncol(module_mat), function(i){
  model_formula <- reformulate('lineage', response=colnames(module_mat)[i])
  model_data <- cbind(module_mat[,i,drop=F], model_frame)
  return(fit_glm(formula=model_formula, data=model_data))
})
names(module_lin_fit_list) <- colnames(module_mat)
lin_model_coefs <- map_dfr(module_lin_fit_list, function(x) x$coefs, .id='module')
lin_model_gof <- map_dfr(module_lin_fit_list, function(x) x$gof, .id='module')

module_pt_fit_list <- pblapply(1:ncol(module_mat), function(i){
  model_formula <- reformulate('pt', response=colnames(module_mat)[i])
  model_data <- cbind(module_mat[,i,drop=F], model_frame)
  return(fit_glm(formula=model_formula, data=model_data))
})
names(module_pt_fit_list) <- colnames(module_mat)
pt_model_coefs <- map_dfr(module_pt_fit_list, function(x) x$coefs, .id='module')
pt_model_gof <- map_dfr(module_pt_fit_list, function(x) x$gof, .id='module')


comb_model_df <- inner_join(comb_model_coefs, comb_model_gof)
lin_model_df <- inner_join(lin_model_coefs, lin_model_gof)
pt_model_df <- inner_join(pt_model_coefs, pt_model_gof)


gof <- bind_rows('lin_pt'=comb_model_gof, 'lin'=lin_model_gof, 'pt'=pt_model_gof, .id='vars') %>% 
  mutate(module=str_replace_all(module, '_', '-'))
coefs <- bind_rows('lin_pt'=comb_model_coefs, 'lin'=lin_model_coefs, 'pt'=pt_model_coefs, .id='vars') %>% 
  mutate(module=str_replace_all(module, '_', '-'))
coefs <- coefs %>% mutate(term=str_remove_all(term, 'lineage'))


# gof %>% write_tsv('data/RNA_ATAC/grn/module_lineage_pt_r2.tsv')
# coefs %>% write_tsv('data/RNA_ATAC/grn/module_lineage_pt_coefs.tsv')

mean_lineage_score <- aggregate_matrix(module_mat, model_frame$lineage)
module_max_lineage <- rownames(mean_lineage_score)[apply(mean_lineage_score, 2, which.max)]
names(module_max_lineage) <- colnames(mean_lineage_score)
module_max_lineage <- enframe(module_max_lineage, 'module', 'module_lineage') %>% 
  mutate(module=str_replace_all(module, '_', '-'))


# gof <- read_tsv('data/RNA_ATAC/grn/module_lineage_pt_r2.tsv')
# coefs <- read_tsv('data/RNA_ATAC/grn/module_lineage_pt_coefs.tsv')



gof_meta <- gof %>% 
  inner_join(gene_meta, by=c('module'='gene')) %>% 
  inner_join(module_max_lineage) %>% 
  mutate(vars=factor(vars, levels=c('lin', 'pt', 'lin_pt')))

ggplot(gof_meta, aes(vars, dsq, color=module_lineage)) +
  geom_quasirandom(alpha=0.5, size=0.2) +
  geom_boxplot(color='black', alpha=0.5, outlier.shape = NA, width=0.4, size=0.2) +
  scale_color_manual(values=lineage_colors) +
  scale_x_discrete(labels=c('lin'='Branch', 'pt'='Pseudotime', 'lin_pt'='Branch +\nPseudotime')) +
  rotate_x_text(40) +
  scale_axis_rangeframe() + theme_rangeframe() + no_legend() + article_text() +
  labs(y=bquote(R^2)) +
  theme(
    axis.title.x = element_blank()
  )

ggsave('plots/paper/fig2/module_r2_lin_pt_boxplot.pdf', width=3.3, height=3.8, unit='cm')


tf_meta <- motif2tf %>% 
  distinct(tf, exp_evidence)
gof_wide <- gof_meta %>% 
  pivot_wider(names_from='vars', values_from='dsq') %>% 
  inner_join(tf_meta, by=c('module'='tf'))

ggplot(gof_wide, aes(lin, pt, fill=module_lineage)) +
  geom_point(data=filter(gof_wide, exp_evidence!='direct'), shape=21, stroke=0.1, size=0.8, fill='grey', alpha=0.5) +
  geom_point(data=filter(gof_wide, exp_evidence=='direct'), shape=21, stroke=0.1, size=0.8) +
  scale_fill_manual(values=lineage_colors) +
  scale_axis_rangeframe() + theme_rangeframe() + no_legend() + article_text() +
  labs(x='Variance explained by branch', y='Variance explained by pseudotime') 

ggsave('plots/paper/fig2/module_r2_lin_pt_scatter.pdf', width=5.2, height=3.8, unit='cm')


ggplot(gof_wide, aes(lin, pt, fill=module_lineage, color=module_lineage, label=module)) +
  geom_text(data=filter(gof_wide, exp_evidence=='direct')) +
  geom_point(data=filter(gof_wide, exp_evidence=='direct'), shape=21, stroke=0.1, size=0.8) +
  scale_fill_manual(values=lineage_colors) +
  scale_color_manual(values=lineage_colors) +
  scale_axis_rangeframe() + theme_rangeframe() + no_legend() + article_text() +
  labs(x='Variance explained by branch', y='Variance explained by pseudotime') 




#### Feature plots of module expression and activity #####
goi <- c('SCRT1', 'NR1D1', 'EMX1', 'IRX5', 'ELK3', 'SP9', 'CCDC88A', 'NEUROD6')

rnatac@active.assay <- 'RNA'
p1 <- feature_plot(rnatac, features = goi, order=T, ncol=length(goi)) &
  scale_color_gradientn(colors=gyylgnbu())
rnatac@active.assay <- 'module_score'
p2 <- feature_plot(rnatac, features = goi, order=T, ncol=length(goi)) &
  scale_color_gradientn(colors=grad(pals::ocean.tempo, 0.5))

p1 / p2


#### Compare motifs ####
gne <- 'HES4'
motif2tf %>% filter(tf==gne)
moi <- motif2tf %>% filter(tf==gne) %>% pull(motif) 
MotifPlot(rnatac, moi[2], assay='peaks')


#### Frac unvalidated ####

motif2tf <- motif2tf %>% 
  mutate(exp_evidence=case_when(
    (origin=='JASPAR2020' & collection != 'UNVALIDATED') |
      (origin=='CIS-BP' & status=='D') ~ 'direct',
    (origin=='JASPAR2020' & collection == 'UNVALIDATED') |
      (origin=='CIS-BP' & status=='I') ~ 'indirect',
    origin == 'SEQ_SIMILARITY' ~ 'inferred'
  )) 


plot_df <- motif2tf %>% 
  distinct(origin, exp_evidence, tf) %>% 
  mutate(origin=factor(origin, levels=c('JASPAR2020', 'CIS-BP', 'SEQ_SIMILARITY')))
ggplot(plot_df, aes(origin, fill=exp_evidence)) +
  geom_bar(color='black', size=0.2) +
  theme_rangeframe() + scale_axis_rangeframe() + 
  article_text() +
  scale_y_continuous(breaks=c(0, 500, 1000), limits=c(0,1000)) +
  rotate_x_text(40) +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y='Count')





#### Same for neg modules ####

module_mat <- t(rnatac[['module_neg_score']]@data[, cells_use])

# module_mat <- t(rnatac[['RNA']]@data[rownames(rnatac[['module_score']]), cells_use])
colnames(module_mat) <- str_replace(colnames(module_mat), '-', '_')

module_fit_list <- pblapply(1:ncol(module_mat), function(i){
  model_formula <- reformulate(colnames(model_frame), response=colnames(module_mat)[i])
  model_data <- cbind(module_mat[,i,drop=F], model_frame)
  return(fit_glm(formula=model_formula, data=model_data))
})
names(module_fit_list) <- colnames(module_mat)
comb_model_coefs <- map_dfr(module_fit_list, function(x) x$coefs, .id='module')
comb_model_gof <- map_dfr(module_fit_list, function(x) x$gof, .id='module')

module_lin_fit_list <- pblapply(1:ncol(module_mat), function(i){
  model_formula <- reformulate('lineage', response=colnames(module_mat)[i])
  model_data <- cbind(module_mat[,i,drop=F], model_frame)
  return(fit_glm(formula=model_formula, data=model_data))
})
names(module_lin_fit_list) <- colnames(module_mat)
lin_model_coefs <- map_dfr(module_lin_fit_list, function(x) x$coefs, .id='module')
lin_model_gof <- map_dfr(module_lin_fit_list, function(x) x$gof, .id='module')

module_pt_fit_list <- pblapply(1:ncol(module_mat), function(i){
  model_formula <- reformulate('pt', response=colnames(module_mat)[i])
  model_data <- cbind(module_mat[,i,drop=F], model_frame)
  return(fit_glm(formula=model_formula, data=model_data))
})
names(module_pt_fit_list) <- colnames(module_mat)
pt_model_coefs <- map_dfr(module_pt_fit_list, function(x) x$coefs, .id='module')
pt_model_gof <- map_dfr(module_pt_fit_list, function(x) x$gof, .id='module')


comb_model_df <- inner_join(comb_model_coefs, comb_model_gof)
lin_model_df <- inner_join(lin_model_coefs, lin_model_gof)
pt_model_df <- inner_join(pt_model_coefs, pt_model_gof)


gof <- bind_rows('lin_pt'=comb_model_gof, 'lin'=lin_model_gof, 'pt'=pt_model_gof, .id='vars') %>% 
  mutate(module=str_replace_all(module, '_', '-'))
coefs <- bind_rows('lin_pt'=comb_model_coefs, 'lin'=lin_model_coefs, 'pt'=pt_model_coefs, .id='vars') %>% 
  mutate(module=str_replace_all(module, '_', '-'))
coefs <- coefs %>% mutate(term=str_remove_all(term, 'lineage'))


# gof %>% write_tsv('data/RNA_ATAC/grn/module_neg_lineage_pt_r2.tsv')
# coefs %>% write_tsv('data/RNA_ATAC/grn/module_neg_lineage_pt_coefs.tsv')





#### Graph viz of branch specific genes ####

specific_tfs <- gene_meta %>% 
  filter(gene%in%names(module_pos)) %>% 
  filter(lineage%in%c('DF', 'VF')) %>% 
  mutate(lineage=fct_recode(lineage, 'ctx'='DF', 'ge'='VF'))

de_modules <- module_score_second_de %>% 
  inner_join(gene_meta, by=c('feature'='gene')) %>%
  # inner_join(specific_tfs, by=c('feature'='gene')) %>% 
  # filter(group==lineage) %>% 
  filter(fc>0, feature%in%tf_graph_meta$name, padj<0.05) %>% 
  group_by(group) %>% 
  top_n(15, fc)

plot_graph <- tf_graph %N>% 
  select(name) %>% 
  inner_join(de_modules, by=c('name'='feature')) %>% 
  mutate(pt_plot=ifelse(group=='ge', -mean_ranked_pt, mean_ranked_pt)) %E>% 
  mutate(forward=.N()$pt_plot[from] < .N()$pt_plot[to]) %>% 
  reroute(
    from=case_when(
      forward & (estimate < 0) ~ to,
      !forward & (estimate > 0) ~ to,
      T ~ from
    ),
    to=case_when(
      forward & (estimate < 0) ~ from,
      !forward & (estimate > 0) ~ from,
      T ~ to
    )
  ) %>% 
  mutate(padj=pmax(padj, 1e-10)) 

ggraph(plot_graph, layout='hive', axis=group, sort.by=mean_ranked_pt) +
  geom_edge_arc(aes(color=sign(estimate)*-log10(padj)), width=0.2) +
  geom_node_text(aes(label=name, color=mean_ranked_pt), size=5/ggplot2::.pt) +
  scale_edge_color_gradientn(colors=pals::ocean.curl(100), limits=c(-10,10)) +
  scale_color_viridis(option='magma', direction=-1, end=0.8) +
  theme_void() + no_legend()


ggsave('plots/paper/fig2/dorsal_ventral_tfs_full2_arc.pdf', width=2.8, height=5.5, unit='cm')  


ggraph(plot_graph, layout='hive', axis=lineage, sort.by=mean_ranked_pt) +
  geom_edge_arc(aes(color=sign(estimate)*-log10(padj)), width=0.2) +
  geom_node_text(aes(label=name), size=5/ggplot2::.pt) +
  scale_edge_color_gradientn(colors=pals::ocean.curl(100), limits=c(-10,10)) +
  theme_void() + no_legend()
ggsave('plots/paper/fig2/dorsal_ventral_tfs_full_arc.pdf', width=2.8, height=5.5, unit='cm')  


ggraph(plot_graph, layout='hive', axis=lineage, sort.by=mean_ranked_pt) +
  geom_axis_hive(label=F, size=0.1) +
  geom_edge_arc(aes(filter=estimate>0), color='#882160', width=0.1) +
  geom_node_point(size=0.01) +
  scale_edge_color_gradientn(colors=pals::ocean.balance(100), limits=c(-10,10)) +
  theme_void()
ggsave('plots/paper/fig2/dorsal_ventral_tfs_pos_arc.pdf', width=1, height=1.7, unit='cm')  



ggraph(plot_graph, layout='hive', axis=lineage, sort.by=mean_ranked_pt) +
  geom_axis_hive(label=F, size=0.1) +
  geom_edge_arc(aes(filter=estimate<0), color='#20877B', width=0.1) +
  geom_node_point(size=0.01) +
  scale_color_viridis(option='magma', direction=-1) +
  scale_edge_color_gradientn(colors=pals::ocean.balance(100), limits=c(-10,10)) +
  theme_void()
ggsave('plots/paper/fig2/dorsal_ventral_tfs_neg_arc.pdf', width=1, height=1.7, unit='cm')  


ggraph(plot_graph, layout='hive', axis=lineage, sort.by=mean_ranked_pt) +
  geom_axis_hive(label=F, size=0.1) +
  geom_edge_arc(aes(filter=estimate<0), color='white') +
  geom_node_point(aes(color=mean_ranked_pt), size=0.02) +
  scale_color_viridis(option='magma', direction=-1) +
  theme_void() +
  no_legend()
ggsave('plots/paper/fig2/dorsal_ventral_tfs_pt_arc.pdf', width=1, height=1.7, unit='cm')  


ggraph(plot_graph, layout='hive', axis=lineage, sort.by=mean_ranked_pt) +
  geom_axis_hive(label=F, size=0.1) +
  geom_edge_arc(aes(filter=estimate<0), color='white') +
  geom_node_point(aes(color=mean_ranked_pt), size=0.02) +
  scale_color_viridis(option='magma', direction=-1) +
  theme_void() +
  no_legend()
ggsave('plots/paper/fig2/dorsal_ventral_tfs_pt_arc.pdf', width=1, height=1.7, unit='cm')  


ggraph(plot_graph, layout='hive', axis=lineage, sort.by=mean_ranked_pt) +
  geom_axis_hive(label=F, size=0.1) +
  geom_edge_arc(aes(filter=estimate<0), color='white') +
  geom_node_point(aes(color=lineage), size=0.02) +
  scale_color_manual(values=lineage_colors) +
  theme_void() +
  no_legend()
ggsave('plots/paper/fig2/dorsal_ventral_tfs_lin_arc.pdf', width=1, height=1.7, unit='cm')  


print_scale(pals::ocean.curl(100))
ggsave('plots/paper/fig2/ocean_curl.pdf') 


ggraph(plot_graph, layout='hive', axis=lineage, sort.by=mean_ranked_pt) +
  geom_axis_hive(label=F, size=0.1) +
  geom_edge_arc(aes(filter=estimate<0, alpha=stat(index)), color='black') +
  geom_node_point(aes(color=mean_ranked_pt), size=0.02) +
  scale_color_viridis(option='magma', direction=-1) +
  theme_void() +
  no_legend()



#### Pull out interaction ####

ctx_modules <- module_score_second_de %>% 
  inner_join(specific_tfs, by=c('feature'='gene')) %>% 
  filter(fc>0, feature%in%tf_graph_meta$name, group==lineage, lineage=='ctx') %>% 
  top_n(15, fc)

ctx_graph <- tf_graph %N>% 
  select(name) %>% 
  inner_join(ctx_modules, by=c('name'='feature')) %E>% 
  filter(estimate>0)

ggraph(ctx_graph) +
  geom_edge_arc(aes(filter=estimate<0, alpha=stat(index)), color='black') +
  geom_node_point(aes(color=mean_ranked_pt), size=2) +
  geom_node_text(aes(label=name)) +
  scale_color_viridis(option='magma', direction=-1) +
  theme_void() +
  no_legend()


ctx_edge_df <- ctx_graph %E>% as_tibble()
from_tf <- 'NEUROD1'
to_tf <- 'EMX1'

edge_use <- ctx_edge_df %>% filter(to_node==to_tf)



modules %>% filter(tf=='GLI3', target%in%de_modules$feature) %>% View


#### Plot tracks for binding peak ####
atac_cov <- subset(atac, lineage%in%c('nt', 'telencephalon', 'early', 'ge', 'ctx'))
atac_cov$lineage <- factor(atac_cov$lineage, levels=c('early', 'nt', 'telencephalon', 'ge', 'ctx'))

to_tf <- 'DLX5'

annot <- gene_annot[gene_annot$gene_name==to_tf]
strand(annot) <- '-'
annot_ex <- Extend(annot, upstream = 100000)
region <- GRanges(seqnames = as.character(seqnames(annot_ex))[[1]], ranges = IRanges(start = min(start(annot_ex)), end = max(end(annot_ex))))
# region <- resize(region, width = 30000, fix='start')


binding_regions <- modules %>% filter(target==to_tf) %>% pull(peaks) %>% str_split(';') %>% unlist() %>% StringToGRanges(sep=c('_', '_'))

binding_regions <- modules %>% filter(target==to_tf) %>% 
  separate_rows(peaks, sep=';') %>% 
  separate(peaks, sep='_', into=c('seqnames', 'start', 'end'), convert=T) %>% 
  filter(end>start(region), start<end(region)) %>% 
  group_by(start, end) %>% 
  mutate(padj_peak=min(padj)) %>% 
  mutate(padj_peak=pmax(padj, 1e-10)) %>% 
  mutate(pos_ratio=sum(sign(estimate)==1)/n()) %>% 
  ungroup()


ggplot(binding_regions) +
  geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
  theme_void() +
  scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
  ggnewscale::new_scale_color() +
  geom_text_repel(
    aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
    size=4, max.overlaps = 100, nudge_y=-100, angle=90, force=5
    
  ) +
  scale_color_manual(values=c('#20877B', '#882160')) +
  xlim(c(start(region), end(region))) + no_legend()


binding_regions_label <- binding_regions %>% 
  filter(padj_peak<1e-3) %>% 
  group_by(start, end) %>% 
  top_n(2, estimate)

ggplot(binding_regions_label) +
  geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
  theme_void() +
  scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
  ggnewscale::new_scale_color() +
  geom_text_repel(
    aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
    size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
      
  ) +
  scale_color_manual(values=c('#20877B', '#882160')) +
  xlim(c(start(region), end(region))) + no_legend()

ggsave('plots/paper/fig2/emx1_reg_sites_labels.pdf', width=8, unit='cm')



p_mot <- ggplot(binding_regions) +
  geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj)), size=4) +
  theme_void() +
  scale_color_gradientn(colors=greys()) +
  xlim(c(start(region), end(region))) +
  no_legend()

p_cov <- CoveragePlot(
  atac_cov, 
  region = region, 
  group.by = 'lineage', 
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=lineage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void() 

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()

p_cov / p_an / p_pe / p_mot + plot_layout(heights=c(10,1,1,1)) & theme(text = element_text(size=5))
ggsave('plots/paper/fig2/emx1_reg_sites_tracks.pdf', width=4.4, height=3.7, unit='cm')

PeakPlot(
  atac_cov,
  region = region
) + scale_x_continuous(breaks=seq(1,1000000000, 1000))
ggsave('plots/paper/fig2/emx1_peaks.pdf', unit='cm', height=3)






#### Tree graph plots with expression and motif ####
cluster_counts <- as.numeric(table(rnatac$RNA_snn_res.20))
cluster_perc_expr <- (GetAssayData(rnatac, assay='RNA', slot='data') > 0) %>% t() %>% 
  aggregate.Matrix(groupings = rnatac$RNA_snn_res.20) %>% {./cluster_counts}

cluster_expr <- cluster_summaries$RNA
cluster_score <- cluster_summaries$module_score
cluster_peaks <- cluster_summaries$module_peaks_chromvar


pick_genes <- c('NEUROD1', 'EMX1')

gene_expr_df <- cluster_expr[, pick_genes] %>% 
  as_tibble(rownames='name')

perc_expr_df <- cluster_perc_expr[, pick_genes] %>% 
  {colnames(.) <- paste0(colnames(.), '_perc');.} %>% 
  as_tibble(rownames='name') 

plot_tree <- lineage_graph %N>% 
  left_join(gene_expr_df) %>% 
  left_join(perc_expr_df)

plots <- map(pick_genes, function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.2) +
    geom_node_point(aes_string(fill=gene, size=paste0(gene, '_perc')), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.5, 3)) +
    scale_fill_gradientn(colors=gyylgnbu()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend() + ggtitle(gene)
})
p_ex <- wrap_plots(plots)



gene_act_df <- cluster_peaks[, pick_genes] %>% 
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>% 
  left_join(gene_act_df)

plots <- map(pick_genes, function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.2) +
    geom_node_point(aes_string(fill=gene, size=gene), size=2, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=gyylorbr()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend() + ggtitle(gene)
})
p_peaks <- wrap_plots(plots)


gene_act_df <- cluster_score[, pick_genes] %>% 
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>% 
  left_join(gene_act_df)

plots <- map(pick_genes, function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.2) +
    geom_node_point(aes_string(fill=gene, size=gene), size=2, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=pals::ocean.tempo(100)) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend() + ggtitle(gene)
})
p_score <- wrap_plots(plots)


p_ex / p_peaks / p_score
ggsave('plots/paper/fig2/neurod6_emx1_expr_peaks_score.png', width=5, height=8)



#### Celltype/lineage-specific GRNs (Analysis from Zhisong) ####

#### Plot graphs ####
sub_grns <- read_rds('/links/groups/treutlein/USERS/zhisong_he/Work/brain_organoid_cropseq/analysis/celltype_GRN/res.lineage_subGRN_dorsal-ventral-nt_comparison_tbl_graph.rds')
sub_grns <- read_rds('/links/groups/treutlein/USERS/zhisong_he/Work/brain_organoid_cropseq/analysis/celltype_GRN/')
edge_enrichment <- read_rds('/links/groups/treutlein/USERS/zhisong_he/Work/brain_organoid_cropseq/analysis/celltype_GRN/res.enrichment_peaks_lineage.rds') %>% 
    as_tibble()

edge_enrichment$region <- as.character(names(grnreg2peaks)[match(edge_enrichment$feature, grnreg2peaks)])


edge_df <- sub_grns$nt %E>% as_tibble()
edge_regions <- edge_df$peak %>% str_replace_all('_', '-')

nt_top_regions <- edge_enrichment %>% 
    filter(group=='nt', region%in%edge_regions, oddsratio>1, padj<0.05) %>% 
    arrange(desc(oddsratio)) %>% filter(row_number()<=100) %>%
    pull(region) %>% str_replace_all('-', '_')

nt_graph <- sub_grns$nt %E>% 
    mutate(is_top=peak%in%nt_top_regions) %>% 
    mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %N>% 
    mutate(prop_nt=num_targets_nt_only/num_targets_global)

ggraph(nt_graph, x=UMAP1, y=UMAP2) +
    geom_edge_diagonal(aes(filter=!is_top), color='gray', width=0.1, alpha=0.1) +
    geom_edge_diagonal(aes(filter=is_top), color='black', width=0.1) +
    geom_node_point(aes(fill=prop_nt_only>0.15, size=num_targets_nt, alpha=prop_nt>0.15), shape=21, stroke=0.5) +
    # geom_node_text(aes(filter=prop_nt_only>0.15, label=name)) +
    scale_alpha_discrete(range=c(0.2,1)) +
    scale_radius(range = c(1,8), trans = 'sqrt') +
    # scale_fill_gradientn(colours = colorRampPalette(c('#ffffff', '#78d5ed', '#045369'))(50), trans = "log1p", na.value = 'white') +
    scale_fill_gradientn(colours = colorRampPalette(c('#78d5ed', '#045369'))(50), trans = "log1p", na.value = 'white') +
    theme_void() + 
    no_legend()
ggsave('plots/paper/fig2/nt_sub_grn_umap.png', width=4, height=4)


ggraph(sub_grns$nt, x=UMAP1, y=UMAP2) +
    geom_node_point(aes(fill=prop_nt_only * 100, size=num_targets_nt), shape=21, stroke=0.5) +
    geom_node_text(aes(label=name))

ggsave('plots/paper/fig2/nt_sub_grn_umap_labels.png', width=20, height=20)



print_scale(colorRampPalette(c('#ffffff', '#78d5ed', '#045369'))(50))
ggsave('plots/paper/fig2/nt_grad.pdf', width=4, height=4)


edge_df <- sub_grns$ctx %E>% as_tibble()
edge_regions <- edge_df$peak %>% str_replace_all('_', '-')

ctx_top_regions <- edge_enrichment %>% 
    filter(group=='ctx', region%in%edge_regions, oddsratio>1, padj<0.05) %>% 
    arrange(desc(oddsratio)) %>% filter(row_number()<=100) %>%
    pull(region) %>% str_replace_all('-', '_')

ctx_graph <- sub_grns$ctx %E>% 
    mutate(is_top=peak%in%ctx_top_regions) %N>% 
    mutate(prop_dorsal=num_targets_dorsal_only/num_targets_global)

ggraph(ctx_graph, x=UMAP1, y=UMAP2) +
    geom_edge_diagonal(aes(filter=!is_top), color='gray', width=0.1, alpha=0.1) +
    geom_edge_diagonal(aes(filter=is_top), color='black', width=0.1) +
    geom_node_point(aes(fill=prop_dorsal>0.15, size=num_targets_dorsal, alpha=prop_dorsal>0.15), shape=21, stroke=0.5) +
    scale_alpha_discrete(range=c(0.2,1)) +
    scale_radius(range = c(1,8), trans = 'sqrt') +
    # scale_fill_gradientn(colours = colorRampPalette(c('#ffffff', '#e65243', '#803028'))(50), trans = "log1p") +
    scale_fill_gradientn(colours = colorRampPalette(c('#e65243', '#803028'))(50), trans = "log1p") +
    theme_void() + 
    no_legend()

ggsave('plots/paper/fig2/ctx_sub_grn_umap.png', width=4, height=4)


print_scale(colorRampPalette(c('#ffffff', '#e65243', '#803028'))(50))
ggsave('plots/paper/fig2/ctx_grad.pdf', width=4, height=4)


edge_df <- sub_grns$ge %E>% as_tibble()
edge_regions <- edge_df$peak %>% str_replace_all('_', '-')

ge_top_regions <- edge_enrichment %>% 
    filter(group=='ge', region%in%edge_regions, oddsratio>1, padj<0.05) %>% 
    arrange(desc(oddsratio)) %>% filter(row_number()<=100) %>%
    pull(region) %>% str_replace_all('-', '_')

ge_graph <- sub_grns$ge %E>% 
    mutate(is_top=peak%in%ge_top_regions) %N>% 
    mutate(prop_ventral=num_targets_ventral_only/num_targets_global)


ggraph(ge_graph, x=UMAP1, y=UMAP2) +
    geom_edge_diagonal(aes(filter=!is_top), color='gray', width=0.1, alpha=0.1) +
    geom_edge_diagonal(aes(filter=is_top), color='black', width=0.1) +
    geom_node_point(aes(fill=prop_ventral>0.1, size=num_targets_ventral, alpha=prop_ventral>0.1), shape=21, stroke=0.5) +
    # geom_node_text(aes(filter=prop_ventral_only>0.15, label=name)) +
    scale_alpha_discrete(range=c(0.2,1)) +
    scale_radius(range = c(1,8), trans = 'sqrt') +
    # scale_fill_gradientn(colours = colorRampPalette(c('#ffffff', '#cc8ae6', '#4a1161'))(50), trans = "log1p") +
    scale_fill_gradientn(colours = colorRampPalette(c('#cc8ae6', '#4a1161'))(50), trans = "log1p") +
    theme_void() + 
    no_legend()

ggsave('plots/paper/fig2/ge_sub_grn_umap.png', width=4, height=4)

print_scale(colorRampPalette(c('#ffffff', '#cc8ae6', '#4a1161'))(50))
ggsave('plots/paper/fig2/ge_grad.pdf', width=4, height=4)


#### Summary stats ####
tf_exp <- motif2tf %>% 
    filter(exp_evidence=='direct') %>% 
    pull(tf) %>% unique()


plot_df <- ge_graph %N>% as_tibble() %>% 
    mutate(num_targets_other=num_targets_global-num_targets_dorsal_only-num_targets_ventral_only-num_targets_nt_only) %>% 
    pivot_longer(c(num_targets_dorsal_only, num_targets_other, num_targets_nt_only, num_targets_ventral_only), names_to='lineage', values_to='targets') %>% 
    arrange(num_targets_global) %>% mutate(name=factor(name, levels=unique(.$name))) %>% 
    filter(pmax(prop_nt_only, prop_dorsal_only, prop_ventral_only)>0.20) %>% 
    mutate(lineage=factor(fct_recode(lineage, 'ctx'='num_targets_dorsal_only', 'ge'='num_targets_ventral_only', 'nt'='num_targets_nt_only', 'early'='num_targets_other'), levels=c('nt', 'ctx', 'ge', 'early'))) %>% 
    group_by(name) %>% mutate(lin_max=factor(which.max(targets[lineage!='early'])), levels=c('nt', 'ctx', 'ge', 'early')) %>% 
    filter(name%in%tf_exp) 


ggplot(plot_df, aes(name, targets, fill=lineage)) +
    geom_bar(stat='identity') +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_fill_manual(values=lineage_colors) +
    coord_flip() +
    facet_grid(lin_max~., space='free', scales='free') +
    article_text() +
    no_legend() +
    labs(x='Transcription factor', y='# targets') +
    no_y_text() +
    no_margin() +
    theme(
        axis.line.y = element_blank(),
        panel.spacing = unit(0, 'lines')
    )

ggsave('plots/paper/fig2/sub_grn_tf_bar.pdf', width=2.6, height=4.8, unit='cm')


ggplot(plot_df, aes(name, targets, fill=lineage)) +
    geom_bar(stat='identity') +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_fill_manual(values=lineage_colors) +
    coord_flip() +
    facet_grid(lin_max~., space='free', scales='free') +
    article_text() +
    no_legend() +
    labs(x='Transcription factor', y='# targets') +
    no_margin() +
    theme(
        axis.line.y = element_blank(),
        panel.spacing = unit(0, 'lines')
    )






#### Plot edges as tracks ####

#### NT ####
atac_cov <- subset(atac, lineage%in%c('nt', 'ctx', 'ge'))
atac_cov$lineage <- factor(atac_cov$lineage, levels=c('nt', 'ge', 'ctx'))

module_peaks <- modules %>% 
    separate_rows(peaks, sep=';') %>% pull(peaks) %>% str_replace_all('_', '-')

enriched_edges <- edge_enrichment %>% filter(group=='nt', !is.na(region), padj<0.01, region%in%module_peaks) 
enriched_regions <- enriched_edges$region %>% unique() %>% str_replace_all('-', '_')

enriched_module_peaks <- modules %>% 
    separate_rows(peaks, sep=';') %>% 
    filter(peaks%in%enriched_regions)

#### Individual TF-TF edge ####
to_tf <- 'ARID5B'
nt_edges <- nt_graph %E>% as_tibble() %>% filter(is_top) %>% select(from_node, to_node, term, peak) %>% 
    filter(to_node%in%to_tf)

annot <- gene_annot[gene_annot$gene_name==to_tf]
strand(annot) <- '-'
annot_ex <- Extend(annot, upstream = 50000)
region <- GRanges(seqnames = as.character(seqnames(annot_ex))[[1]], ranges = IRanges(start = min(start(annot_ex)), end = max(end(annot_ex))))
# region <- resize(region, width = 30000, fix='start')
binding_regions <- nt_edges %>% 
    separate(peak, sep='_', into=c('seqnames', 'start', 'end'), convert=T) %>% 
    filter(end>start(region), start<end(region))


p_mot <- ggplot(binding_regions) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys()) +
    xlim(c(start(region), end(region))) +
    no_legend()

p_cov <- CoveragePlot(
    atac_cov, 
    region = region, 
    group.by = 'lineage', 
    annotation = FALSE,
    peaks = FALSE
) & scale_fill_manual(values=lineage_colors) &
    theme_void() & no_legend() & theme(strip.text = element_blank())

p_an <- AnnotationPlot(
    atac_cov,
    region = region
) + theme_void() 

p_pe <- PeakPlot(
    atac_cov,
    region = region
) + theme_void() + no_legend()

p_cov / p_an / p_pe / p_mot + plot_layout(heights=c(10,1,1,1)) & theme(text = element_text(size=5))


#### ISL1 ####


to_tf <- 'PAX6'
annot <- gene_annot[gene_annot$gene_name==to_tf]
strand <- mode(as.character(strand(annot)))
# strand <- '+'
region <- GRanges(
    seqnames = as.character(seqnames(annot))[[1]], 
    ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
    strand = strand
)
# region <- resize(region, width=0, fix='start')
region <- Extend(region, upstream = 50000, downstream=10)

binding_regions <- modules %>% filter(target==to_tf) %>% pull(peaks) %>% str_split(';') %>% unlist() %>% StringToGRanges(sep=c('_', '_'))

binding_regions <- modules %>% filter(target==to_tf) %>% 
    separate_rows(peaks, sep=';') %>% 
    filter(peaks%in%enriched_regions) %>% 
    separate(peaks, sep='_', into=c('seqnames', 'start', 'end'), convert=T) %>% 
    filter(end>start(region), start<end(region)) %>% 
    group_by(start, end) %>% 
    mutate(padj_peak=min(padj)) %>% 
    mutate(padj_peak=pmax(padj, 1e-10)) %>% 
    mutate(pos_ratio=sum(sign(estimate)==1)/n()) %>% 
    ungroup()


ggplot(binding_regions) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=4, max.overlaps = 100, nudge_y=-100, angle=90, force=5
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


binding_regions_label <- binding_regions %>% 
    # filter(padj_peak<1e-4) %>% 
    group_by(start, end) %>% 
    top_n(3, -padj_peak)

ggplot(binding_regions_label) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


plot_df <- binding_regions %>% 
    group_by(start, end) %>% 
    summarize(padj_peak=min(padj_peak))
p_mot <- ggplot(plot_df) +
    geom_text_repel(
        data = binding_regions_label,
        mapping = aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    ggnewscale::new_scale_color() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    xlim(c(start(region), end(region))) + no_legend()

p_cov <- CoveragePlot(
    atac_cov, 
    region = region, 
    group.by = 'lineage', 
    annotation = FALSE,
    peaks = FALSE
) & scale_fill_manual(values=lineage_colors) &
    theme_void() & no_legend() & theme(strip.text = element_blank())

p_an <- AnnotationPlot(
    atac_cov,
    region = region
) + theme_void() 

p_pe <- PeakPlot(
    atac_cov,
    region = region
) + theme_void() + no_legend()

p_cov / p_an / p_pe / p_mot + plot_layout(heights=c(5,1,1,3)) & theme(text = element_text(size=5))



dim_plot(rnatac, group.by='stage_manual')


#### Look at regions ####


ct_grn <- read.table('/links/groups/treutlein/USERS/zhisong_he/Work/brain_organoid_cropseq/analysis/celltype_GRN/res.GRN_trimming_lineages.tsv') %>% 
    as_tibble()

ct_grn <- read_tsv('data/RNA_ATAC/grn/networks/v5_d100k_cons_cre_noex_ct_specific_grn.tsv')
# ct_grn %>% write_tsv('data/RNA_ATAC/grn/networks/v5_d100k_cons_cre_noex_ct_specific_grn.tsv')

ct_grn <- ct_grn %>% select(tf, target, corr, log_padj, estimate, padj, peak, ctx:telencephalon) %>% 
    filter(tf%in%tf_exp)
plot_df <- ct_grn %>% 
    pivot_longer(ctx:telencephalon) %>% 
    group_by(tf, name) %>% 
    summarize(n_con=sum(value))

ggplot(plot_df, aes(name, log(n_con), fill=name)) +
    geom_boxplot()

plot_df <- ct_grn %>% 
    pivot_longer(ctx:telencephalon) %>% 
    group_by(tf, peak) %>% 
    summarize(
        frac_active=sum(value)/n()
    ) %>% 
    filter(frac_active>0) %>% 
    group_by(tf) %>% 
    summarize(
        frac_active = mean(frac_active),
        n_regions = n()
    )

ggplot(plot_df, aes(log10(n_regions), frac_active, label=tf)) +
    geom_text(size=2)


ggsave('plots/RNA_ATAC/grn/frac_active_tf_scatter.png', width=8, height=5)


#### Get mean expression of TFs in lineages ####

tf_lineage_expr <- aggregate_matrix(t(GetAssayData(rnatac)[unique(ct_grn$tf), ]), groups=rnatac$lineage)
tf_expr_df <- tf_lineage_expr %>% as_tibble(rownames='lineage') %>% 
    pivot_longer(!lineage, names_to='tf', values_to='expr')


plot_df <- ct_grn %>% 
    filter(ctx | ge | nt) %>% 
    group_by(tf) %>% 
    pivot_longer(ctx:telencephalon, names_to='lineage', values_to='is_active') %>% 
    filter(lineage%in%c('ctx', 'ge')) %>% 
    inner_join(tf_expr_df) %>% 
    mutate(tf_act=expr*estimate) %>% 
    # mutate(tf_act = tf_act/sd(tf_act)) %>% 
    filter(is_active, !any(is.infinite(tf_act))) %>% 
    group_by(tf, lineage) %>% 
    summarize(tf_act=mean(tf_act)) %>% 
    group_by(tf) %>% 
    # filter(n()>1) %>%
    mutate(
        max_est=max(abs(tf_act)),
        span=max(tf_act) - min(tf_act),
        crossing=sign(min(tf_act))!=sign(max(tf_act)),
        abs_span=abs(max(tf_act)) + abs(min(tf_act))
    ) %>% 
    ungroup() %>% 
    group_by(tf) %>% 
    mutate(high_lin=lineage[which.max(abs(tf_act[lineage%in%c('ctx', 'ge')]))])


ctx_high_tfs <- plot_df %>% filter(high_lin=='ctx', lineage=='ctx') %>% mutate(max_act=-max(tf_act)) %>% arrange(max_act) %>% pull(tf)
ge_high_tfs <- plot_df %>% filter(high_lin=='ge', lineage=='ge') %>% mutate(max_act=max(tf_act)) %>% arrange(max_act) %>% pull(tf)
# cross_tfs <- plot_df %>% filter(crossing) %>% arrange(high_lin) %>% pull(tf)


tf_order <- unique(c(ctx_high_tfs, ge_high_tfs))
    

ggplot(plot_df, aes(factor(tf, levels=tf_order), tf_act)) +
    geom_line() +
    geom_point(mapping=aes(color=lineage)) +
    geom_text(mapping=aes(label=tf), size=, angle=90) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values=lineage_colors) +
    theme(axis.text.x = element_text(size=5)) +
    rotate_x_text(40)

ggsave('plots/paper/fig2/tf_lineage_diff_names.pdf', width=40, height=6)


ggplot(plot_df, aes(factor(tf, levels=tf_order), tf_act)) +
    geom_line() +
    geom_point(mapping=aes(color=lineage)) +
    geom_text(data=filter(plot_df, tf%in%target_names), mapping=aes(label=tf), size=2, angle=90) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values=lineage_colors) +
    theme(axis.text.x = element_text(size=5)) +
    rotate_x_text(40)

ggsave('plots/RNA_ATAC/grn/tf_lineage_target_names.png', width=30, height=5)

ggplot(plot_df, aes(factor(tf, levels=tf_order), tf_act)) +
    geom_hline(yintercept = 0, color='darkgrey', linetype='dashed', size=0.1) +
    geom_line(size=0.1) +
    geom_point(mapping=aes(color=lineage), size=0.01, shape=16) +
    scale_color_manual(values=lineage_colors) +
    theme(axis.text.x = element_text(size=5)) +
    scale_axis_rangeframe() +
    theme_rangeframe() +
    article_text() +
    no_x_text() +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    ) +
    no_legend() 

ggsave('plots/paper/fig2/tf_lineage_diff.pdf', width=6.4, height=3, unit='cm')















