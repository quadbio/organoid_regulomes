source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')

setwd('~/projects/early/')

library(ggraph)
library(tidytree)
library(fastmatch)
library(biomaRt)

select <- dplyr::select

#### Plots for figure 1 ####
#### Load data ####
rnatac <- read_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2.1_srt.rds')

feature_sets <- read_rds('data/gene_sets/RNA_feature_sets.rds')

gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')
gene_annot <- gene_annot[gene_annot$gene_name %in% feature_sets$grn_features]
peak_ranges <- StringToGRanges(rownames(rnatac@assays$peaks))

atac <- read_rds('data/ATAC/ATAC_all_merged/ATAC_all_merged_v4_srt.rds')


## Select genes
mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl')
all_coding_genes <- getBM(attributes = c( 'hgnc_symbol'), filters = c('biotype'), values = list(biotype='protein_coding'), mart = mart)
good_genes <- all_coding_genes$hgnc_symbol
good_genes <- good_genes[!str_detect(good_genes, '^MT-|^RP|^HIST')]


### Merge lineage info with atac
ps_matches <- read_tsv('data/RNA_ATAC/integration/all_matches.tsv')
ps_matches_atac <- ps_matches %>%
  distinct(pseudocell, ATAC) %>%
  column_to_rownames('ATAC')
ps_order <- ps_matches_atac[colnames(atac), ]
atac_meta <- rnatac@meta.data[ps_order, c('RNA_snn_res.20', 'lineage', 'lineage_coarse')]
atac_meta <- rownames_to_column(atac_meta, 'pseudocell')
rownames(atac_meta) <- colnames(atac)
atac <- AddMetaData(atac, atac_meta)

motif2tf <- read_tsv('data/tf_matches/motif2tf_early_all.tsv') %>%
  filter(origin=='JASPAR2020', collection=='CORE')

lineage_graph <- read_rds('data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')
lineage_tree <- read_rds('data/RNA_ATAC/lineages/cellrank_lineage_tree.rds')
hc_meta <- read_tsv('data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv') %>%
  mutate(cluster=as.character(cluster))
cluster_meta <- as_tibble(lineage_graph)
# cluster_meta %>% write_tsv('data/RNA_ATAC/lineages/cellrank_graph_res20_meta.tsv')

gene_meta <- read_tsv('data/gene_sets/gene_scores.tsv')

#### General plots of timecourse ####
umap_meta <- Reductions(rnatac, slot='umap')@cell.embeddings %>%
  as_tibble(rownames='cell')

meta <- rnatac@meta.data %>%
  as_tibble(rownames='cell') %>%
  inner_join(umap_meta)


ggplot(meta, aes(x=UMAP_1, y=UMAP_2,  color=factor(age))) +
    geom_point(alpha=0.6, size=1) +
    scale_color_manual(values=age_colors) +
    theme_void() +
    no_legend()
ggsave('plots/paper/fig1/age_umap.png', width=7, height=4)


ggplot(meta, aes(x=UMAP_1, y=UMAP_2,  color=factor(line))) +
    geom_point(alpha=0.6, size=1) +
    scale_color_manual(values=line_colors) +
    theme_void() +
    # facet_grid(line~.) +
    no_legend()
ggsave('plots/paper/fig1/line_umap.png', width=7, height=4)


ggplot(meta, aes(x=UMAP_1, y=UMAP_2,  color=factor(line))) +
    geom_point(alpha=0.6, size=1) +
    scale_color_manual(values=line_colors) +
    theme_void() +
    facet_grid(line~.) +
    no_legend()
ggsave('plots/paper/fig1/line_split_umap.png', width=7, height=16)


ggplot(meta, aes(x=UMAP_1, y=UMAP_2,  color=factor(stage))) +
  geom_point(alpha=0.6, size=1) +
  scale_color_manual(values=stage_colors) +
  theme_void() +
  # facet_grid(line~.) +
  no_legend()
ggsave('plots/paper/fig1/stage_umap.png', width=7, height=4)


#### Feature plots ####
add_markers <- c('ZIC1', 'ZIC4', 'NEUROD4', 'NEUROD2', 'GAD2', 'HOPX', 'EMX2', 'EMX1', 'DLX5')
feature_plot(rnatac, features=paste0('geneactivity_', add_markers), order=T) &
  scale_color_gradientn(colors=gypurd(1))


feature_plot(rnatac, features='NFIB', order=T) &
  scale_color_gradientn(colors=gyylgnbu())


plot_markers <- union(all_markers, add_markers)

print_scale(gybupu(1.5))

feature_plot(rnatac, features=paste0('geneactivity_', plot_markers), order=T) &
  scale_color_gradientn(colors=gypurd())
ggsave('plots/paper/fig1/gene_act_all_feature_umap.png', width=12, height=24)


motifs_use <- motif2tf %>% filter(tf%in%plot_markers)
motifs_use <- motif2tf %>% filter(tf=='NFIB')
plot_motifs <- motifs_use$motif %>% unique()

feature_plot(rnatac, features=paste0('chromvar_', plot_motifs), order=T) &
  scale_color_gradientn(colors=gyylorbr())
ggsave('plots/paper/fig1/chromvar_all_feature_umap.png', width=12, height=16)

feature_plot(rnatac, features=plot_markers, order=T) &
  scale_color_gradientn(colors=gyylgnbu())
ggsave('plots/paper/fig1/rna_all_feature_umap.png', width=12, height=24)


print_scale(gypurd())
ggsave('plots/paper/fig1/gypurd.pdf')

print_scale(gyylorbr())
ggsave('plots/paper/fig1/gyylorbr.pdf')

print_scale(gyylgnbu())
ggsave('plots/paper/fig1/gyylgnbu.pdf')


#### Motif plots ####
MotifPlot(rnatac, assay='peaks', motifs='MA1115.1') +
  scale_fill_manual(values=c(green, blue, yellow,red))
MotifPlot(rnatac, assay='peaks', motifs='MA0069.1') +
  scale_fill_manual(values=c(green, blue, yellow,red))
MotifPlot(rnatac, assay='peaks', motifs='MA0612.2') +
  scale_fill_manual(values=c(green, blue, yellow,red))
MotifPlot(rnatac, assay='peaks', motifs='MA1476.1') +
  scale_fill_manual(values=c(green, blue, yellow,red))


#### Cluster summaries ####
cluster_counts <- as.numeric(table(rnatac$RNA_snn_res.20))
cluster_expr <- GetAssayData(rnatac, assay='RNA', slot='data') %>% t() %>%
  aggregate.Matrix(groupings = rnatac$RNA_snn_res.20) %>% {./cluster_counts}

cluster_perc_expr <- (GetAssayData(rnatac, assay='RNA', slot='data') > 0) %>% t() %>%
  aggregate.Matrix(groupings = rnatac$RNA_snn_res.20) %>% {./cluster_counts}

cluster_motifs <- GetAssayData(rnatac, assay='chromvar', slot='data') %>% t() %>%
  aggregate.Matrix(groupings = rnatac$RNA_snn_res.20) %>% {./cluster_counts}

cluster_act <- GetAssayData(rnatac, assay='gene_activity', slot='data') %>% t() %>%
  aggregate.Matrix(groupings = rnatac$RNA_snn_res.20) %>% {./cluster_counts}

cluster_peaks <- GetAssayData(rnatac, assay='peaks', slot='data') %>% t() %>%
  aggregate.Matrix(groupings = rnatac$RNA_snn_res.20) %>% {./cluster_counts}

cluster_perc_acc <- (GetAssayData(rnatac, assay='peaks', slot='data') > 0) %>% t() %>%
  aggregate.Matrix(groupings = rnatac$RNA_snn_res.20) %>% {./cluster_counts}

# cluster_summaries <- list()
cluster_summaries[['RNA']] <- cluster_expr
cluster_summaries[['perc_expr']] <- cluster_perc_expr
cluster_summaries[['chromvar']] <- cluster_motifs
cluster_summaries[['gene_activity']] <- cluster_act
cluster_summaries[['peaks']] <- cluster_peaks
cluster_summaries[['perc_acc']] <- cluster_perc_acc

# cluster_summaries %>% write_rds('data/RNA_ATAC/RNA_ATAC_res20_summaries.rds')

#### DE between stages ####

expr_stage_de <- de(rnatac, 'stage', assay='RNA')
motif_stage_de <- de(rnatac, 'stage', assay='chromvar')
act_stage_de <- de(rnatac, 'stage', assay='gene_activity')

# expr_stage_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_expr_stage.tsv')
# motif_stage_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_motif_stage.tsv')
# act_stage_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_act_stage.tsv')

# expr_stage_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_expr_stage.tsv')
# motif_stage_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_motif_stage.tsv')
# act_stage_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_act_stage.tsv')



#### Chromatin accessibility dynamics over timecourse ####
cluster_summaries <- read_rds('data/RNA_ATAC/RNA_ATAC_res20_summaries.rds')


de_expr <- expr_stage_de %>%
  group_by(feature) %>%
  filter(padj<1e-4, prcex_self>0.3)


de_act <- act_stage_de %>%
  filter(padj<1e-1)

top_stage_act <- act_stage_de %>%
  group_by(feature) %>%
  filter(abs(fc)==max(abs(fc)), padj<0.01) %>%
  group_by(group) %>% filter(feature%in%de_expr$feature) %>%
  arrange(factor(group, levels=c('EB', 'nepi', 'nect', 'orgnaoid'))) %>%
  top_n(100, fc) %>% pull(feature, group)

top_stage_expr <- expr_stage_de %>%
  group_by(feature) %>%
  filter(abs(fc)==max(abs(fc)), padj<0.01) %>%
  group_by(group) %>% filter(feature%in%de_act$feature) %>%
  arrange(factor(group, levels=c('EB', 'nepi', 'nect', 'orgnaoid'))) %>%
  top_n(100, fc) %>% pull(feature, group)


top_de_fests <- intersect(top_stage_act, top_stage_expr)

top_de_act <- cluster_summaries$gene_activity[, top_de_fests]
top_de_expr <- cluster_summaries$RNA[, top_de_fests]

act_df <- top_de_act %>%
  as_tibble(rownames='name') %>%
  filter(name!=25) %>%
  filter(name!=171) %>%
  pivot_longer(!name, names_to='gene', values_to='expr') %>%
  inner_join(gene_meta) %>% select(-lineage) %>% inner_join(cluster_meta) %>%
  group_by(gene) %>% mutate(expr=scale01(expr)) %>% ungroup() %>%
  mutate(gene=factor(gene, levels=rev(unique(.$gene)))) %>%
  # mutate(gene=factor(gene, levels=genes_use)) %>%
  arrange(pseudotime_ranks) %>%
  mutate(
    name=factor(name, levels=unique(.$name)),
    lineage=factor(lineage, levels=c('early', 'nt', 'telencephalon', 'ge', 'ctx'))
  )

expr_df <- top_de_expr %>%
  as_tibble(rownames='name') %>%
  filter(name!=25) %>%
  filter(name!=171) %>%
  pivot_longer(!name, names_to='gene', values_to='expr') %>%
  inner_join(gene_meta) %>% select(-lineage) %>% inner_join(cluster_meta) %>%
  group_by(gene) %>% mutate(expr=scale01(expr)) %>% ungroup() %>%
  mutate(gene=factor(gene, levels=rev(unique(.$gene)))) %>%
  # mutate(gene=factor(gene, levels=genes_use)) %>%
  arrange(pseudotime_ranks) %>%
  mutate(
    name=factor(name, levels=unique(.$name)),
    lineage=factor(lineage, levels=c('early', 'nt', 'telencephalon', 'ge', 'ctx'))
  )


p_act <- ggplot(act_df, aes(name, gene, fill=expr)) +
  geom_tile() +
  facet_grid(scales='free', space='free') +
  # scale_fill_gradientn(colors=pals::ocean.curl(100)) +
  scale_fill_gradientn(colors=pals::brewer.greys(100)) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0.03, b = 0.03, r = 0.01, l = 0.01, unit='cm'),
    axis.text.y = element_text()
  ) +
  ggtitle('Gene activity')
p_act

p_expr <- ggplot(expr_df, aes(name, gene, fill=expr)) +
  geom_tile() +
  facet_grid(scales='free', space='free') +
  # scale_fill_gradientn(colors=pals::ocean.curl(100)) +
  scale_fill_gradientn(colors=pals::brewer.greys(100)) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0.03, b = 0.03, r = 0.01, l = 0.01, unit='cm'),
    axis.text.y = element_text()
  ) +
  ggtitle('Expression')
p_expr


p_expr / p_act



#### Corr between gene act and gene expression ####

rnatac <- FindVariableFeatures(rnatac, nfeatures=4000)
gene_act_features <- rownames(rnatac@assays$gene_activity@data)
var_feats <- rnatac@assays$RNA@meta.features %>%
  as_tibble(rownames='gene') %>%
  filter(!str_detect(gene, '^MT-|^RP|^HIST'), gene%in%gene_act_features) %>%
  top_n(4000, vst.variance.standardized) %>%
  pull(gene) %>% unique()

expr_mat <- t(GetAssayData(rnatac, assay='RNA')[var_feats, ])
act_mat <- t(GetAssayData(rnatac, assay='gene_activity')[var_feats, ])
motif_mat <- t(GetAssayData(rnatac, assay='chromvar'))

gene_annot_var <- gene_annot[gene_annot$gene_name%in%var_feats, ]
peaks_near_gene <- find_peaks_near_genes(peak_ranges, gene_annot_var, distance = 2000, only_tss = FALSE)
peaks2gene <- aggregate.Matrix(t(peaks_near_gene), groupings=colnames(peaks_near_gene))

expr_act_cor <- sparse_cov(expr_mat, act_mat)
expr_act_self_cor <- map_dbl(rownames(expr_act_cor$cor), function(x) expr_act_cor$cor[x,x])
names(expr_act_self_cor) <- rownames(expr_act_cor$cor)
mean_expr <- colMeans(expr_mat)[rownames(expr_act_cor$cor)]
prc_expr <- percent_expression(expr_mat)[rownames(expr_act_cor$cor)]

peak_counts <- map_dbl(seq(nrow(peaks2gene)), function(x){
  length(unique(rownames(gene2peaks)[as.logical(peaks2gene[x, ])]))
})
names(peak_counts) <- rownames(peaks2gene)
peak_counts <- peak_counts[var_feats]


label_points <- c('NEUROD6', 'NFIA', 'POU5F1', 'KRT8', 'BCL11A', 'FABP7', 'STMN4', 'DLX6', 'LIN28A', 'GRID2')


cor_plot_df <- expr_act_self_cor %>%
  enframe('gene', 'cor') %>%
  mutate(mean_expr=mean_expr, prc_expr=prc_expr, n_peaks=peak_counts)



ggplot(filter(cor_plot_df, gene%in%label_points), aes(n_peaks, cor, label=gene)) +
  geom_text()
ggplot(cor_plot_df, aes(n_peaks, cor, label=gene)) +
  geom_text()


ggplot(cor_plot_df, aes(n_peaks, cor, size=prc_expr, label=gene)) +
  # geom_point(shape=21, stroke=0.1) +
  geom_hline(yintercept=0, size=0.25, color='grey') +
  geom_point(shape=21, stroke=0, alpha=0.4, fill='black') +
  geom_point(data=filter(cor_plot_df, gene%in%label_points), shape=21, stroke=0, alpha=0.99, fill='black') +
  scale_size_continuous(range=c(0.2,1.2)) +
  geom_rangeframe(data=tibble(y=c(-0.1, 0.5), x=c(0,400)), mapping=aes(y=y, x=x), inherit.aes=F, size=0.25) +
  scale_y_continuous(limits=c(-0.1, 0.5), breaks=seq(-0.1, 0.5, 0.1), expand=c(0,0.02)) +
  scale_x_continuous(expand=c(0.05,0), breaks=c(0,200,400)) +
  theme(
    panel.border = element_blank(),
    text = element_text(size=5),
    axis.text = element_text(size=5),
    axis.ticks = element_line(size=0.25),
    axis.ticks.length = unit(0.05, 'cm')
  ) +
  no_legend() +
  labs(x='Number of peaks', y='Pearson correlation')

ggsave('plots/paper/fig1/rna_expr_act_cor_scatter.pdf', width=3, height=4.6, unit='cm')

#### DE between lineages ####

expr_stage_de <- de(rnatac, 'stage', assay='RNA')
motif_stage_de <- de(rnatac, 'stage', assay='chromvar')
act_stage_de <- de(rnatac, 'stage', assay='gene_activity')

expr_lin_de <- de(rnatac, 'lineage', assay='RNA')
motif_lin_de <- de(rnatac, 'lineage', assay='chromvar')
act_lin_de <- de(rnatac, 'lineage', assay='gene_activity')
peak_lin_de <- de(rnatac, 'lineage', assay='peaks')

# expr_lin_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_expr_lineage.tsv')
# motif_lin_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_motif_lineage.tsv')
# act_lin_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_act_lineage.tsv')
# peak_lin_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_peaks_lineage.tsv')

# expr_lin_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_expr_lineage.tsv')
# motif_lin_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_motif_lineage.tsv')
# act_lin_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_act_lineage.tsv')
# peak_lin_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_peaks_lineage.tsv')

de_lin_list <- list(
  RNA = expr_lin_de,
  chromvar = motif_lin_de,
  gene_activity = act_lin_de,
  peaks = peak_lin_de
)

de_lin_list %>% write_rds('data/RNA_ATAC/diff_expression/lineage_de.rds')

expr_lin_coarse_de <- de(rnatac, 'lineage_coarse', assay='RNA')
motif_lin_coarse_de <- de(rnatac, 'lineage_coarse', assay='chromvar')
act_lin_coarse_de <- de(rnatac, 'lineage_coarse', assay='gene_activity')
peak_lin_coarse_de <- de(rnatac, 'lineage_coarse', assay='peaks')


# expr_lin_coarse_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_expr_lineage_coarse.tsv')
# motif_lin_coarse_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_motif_lineage_coarse.tsv')
# act_lin_coarse_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_act_lineage_coarse.tsv')
# peak_lin_coarse_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_peaks_lineage_coarse.tsv')

# expr_lin_coarse_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_expr_lineage_coarse.tsv')
# motif_lin_coarse_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_motif_lineage_coarse.tsv')
# act_lin_coarse_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_act_lineage_coarse.tsv')
# peak_lin_coarse_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_peaks_lineage_coarse.tsv')

de_lin_coarse_list <- list(
  RNA = expr_lin_coarse_de,
  chromvar = motif_lin_coarse_de,
  gene_activity = act_lin_coarse_de,
  peaks = peak_lin_coarse_de
)

de_lin_coarse_list %>% write_rds('data/RNA_ATAC/diff_expression/lineage_coarse_de.rds')


#### Heatmaps ####

cluster_summaries <- read_rds('data/RNA_ATAC/RNA_ATAC_res20_summaries.rds')

## Select DE genes
top_expr_stage <- expr_stage_de %>%
  filter(feature%in%good_genes) %>%
  group_by(group) %>% top_n(5, fc)

top_stage_markers <- top_expr_stage %>%
  arrange(factor(group, levels=names(stage_colors)), feature) %>%
  pull(feature) %>% unique()

top_expr_lin <- expr_lin_de %>%
  filter(feature%in%good_genes) %>%
  mutate(expr_ratio=prcex_self/prcex_other) %>%
  filter(!group%in%c('early', 'other'), expr_ratio>2) %>%
  group_by(group) %>% top_n(40, fc)

top_lin_markers <- top_expr_lin %>%
  arrange(factor(group, levels=names(stage_colors))) %>%
  pull(feature) %>% unique()


filter(top_expr_lin, group=='nt')$feature
filter(top_expr_lin, group=='ge')$feature
filter(top_expr_lin, group=='ctx')$feature

## Custom geen list
genes_show <- list(
  # 'EB' = filter(top_expr_stage, group=='EB')$feature,
  'nepi' = filter(top_expr_stage, group=='nepi')$feature,
  'nect' = filter(top_expr_stage, group=='nect')$feature,
  'org' = filter(top_expr_stage, group=='organoid')$feature,
  'nt' = c('OTX1', 'TFPI', 'NRIP3', 'RSPO1', 'WNT7B', 'LHX9', 'LHX1', 'CDH6'),
  'telencephalon' = c('FOXG1', 'DIO3', 'PAX2', 'SIX6', 'DCT'),
  'ge' = c('DLX5', 'DLX6', 'DLX2', 'SOX6', 'GSX2', 'NKX2-1', 'GAD2', 'PCDH19'),
  'ctx' = c('EMX2', 'TBR1', 'NEUROD6', 'NEUROG2', 'NFIB', 'HOPX', 'IFI44L', 'NFIA')
)

gene_groups <- tibble(
  gene = unlist(genes_show, use.names=F),
  gene_group = rep(names(genes_show), map(genes_show, length))
)

# genes_use <- union(top_stage_markers, top_lin_markers) %>% intersect(colnames(cluster_summaries$RNA))
# genes_use <- top_lin_markers %>% intersect(colnames(cluster_summaries$RNA)) %>%
#   union(c('NKX2-1', 'ZIC1', 'ZIC2', 'FOXG1'))

# genes_use <- filter(top_expr_lin, group=='telencephalon')$feature
genes_use <- unique(unlist(genes_show))
gene_expr <- cluster_summaries$RNA[, genes_use]
# gene_expr <- get_features(rnatac, features=top_stage_markers, assay='RNA', slot='data')

ex_df <- gene_expr %>%
  as_tibble(rownames='name') %>%
  filter(name!=171) %>%
  pivot_longer(!name, names_to='gene', values_to='expr') %>%
  inner_join(gene_meta) %>% select(-lineage) %>% inner_join(cluster_meta) %>% inner_join(gene_groups) %>%
  filter(name!=25) %>%
  group_by(gene) %>% mutate(expr=scale01(expr)) %>% ungroup() %>%
  mutate(gene_group=factor(gene_group, levels=names(genes_show))) %>%
  arrange(gene_group, mean_ranked_pt) %>%
  mutate(gene=factor(gene, levels=rev(unique(.$gene)))) %>%
  # mutate(gene=factor(gene, levels=genes_use)) %>%
  arrange(pseudotime_ranks) %>%
  mutate(
    name=factor(name, levels=unique(.$name)),
    lineage=factor(lineage, levels=c('early', 'nt', 'telencephalon', 'ge', 'ctx'))
  )


ggplot(ex_df, aes(name, gene, fill=expr)) +
  geom_tile() +
  facet_grid(gene_group~lineage, scales='free', space='free') +
  # scale_fill_gradientn(colors=pals::ocean.curl(100)) +
  scale_fill_gradientn(colors=greys(1))


## Select DA motifs
top_motif_stage <- motif_stage_de %>%
  inner_join(motif2tf, by=c('feature'='motif')) %>%
  group_by(group) %>% top_n(20, fc)

top_motif_lin <- motif_lin_de %>%
  inner_join(motif2tf, by=c('feature'='motif')) %>%
  group_by(group) %>% top_n(60, fc)


## Custom motif list
motifs_show <- list(
  'nepi' = filter(top_motif_stage, tf%in%c('POU5F1', 'SREBF1', 'NFKB1', 'SP1', 'KLF6'))$feature,
  'nect' = filter(top_motif_stage, tf%in%c('PITX1', 'NR2F1', 'SOX8', 'RXRG', 'PPARD'))$feature,
  'org' = filter(top_motif_stage, tf%in%c('POU4F1', 'RFX2', 'JDP2', 'ATF3', 'LMX1A'))$feature,
  'nt' = filter(top_motif_lin, tf%in%c('OTX2', 'TCF7L2', 'RHOXF1', 'LEF1', 'PAX6'))$feature,
  'telencephalon' = filter(top_motif_lin, tf%in%c('PAX9', 'PAX1', 'NR2F2', 'NR4F1', 'SOX14', 'FOXG1'))$feature,
  'ge' = filter(top_motif_lin, tf%in%c('TCFL7', 'HOXB8', 'HOXA7', 'DLX5', 'MEIS2'))$feature,
  'ctx' = filter(top_motif_lin, tf%in%c('NFIA', 'NEUROD1', 'EMX1', 'TLX1', 'EN1'))$feature
)

motif_groups <- tibble(
  motif = unlist(motifs_show, use.names=F),
  motif_group = rep(names(motifs_show), map(motifs_show, length))
)


# motif_act <- get_features(rnatac, features=top_stage_motifs, assay='chromvar', slot='data')
motifs_use <- unique(filter(top_motif_lin, group=='telencephalon')$feature)
motifs_use <- unique(filter(top_motif_stage, group=='nect')$feature)
motifs_use <- unique(unlist(motifs_show))
motif_act <- cluster_summaries$chromvar[, motifs_use]

motif_df <- motif_act %>%
  as_tibble(rownames='name') %>%
  filter(name!=171, name!=25) %>%
  pivot_longer(!name, names_to='motif', values_to='expr') %>%
  left_join(select(motif2tf, -name)) %>% inner_join(cluster_meta) %>%
  group_by(motif) %>% mutate(expr=scale01(expr)) %>% ungroup() %>%
  inner_join(motif_groups) %>% mutate(motif_group=factor(motif_group, levels=names(motifs_show))) %>%
  arrange(pseudotime_ranks) %>%
  mutate(name=factor(name, levels=unique(.$name))) %>%
  mutate(
    name=factor(name, levels=unique(.$name)),
    lineage=factor(lineage, levels=c('early', 'nt', 'telencephalon', 'ge', 'ctx'))
  )

ggplot(motif_df, aes(name, tf, fill=expr)) +
  geom_tile() +
  facet_grid(motif_group~lineage, scales='free', space='free') +
  scale_fill_gradientn(colors=ylorbr(1.5))


meta <- rnatac@meta.data %>%
  as_tibble(rownames='cell') %>%
  filter(!RNA_snn_res.20%in%c(171,142,25,162)) %>% select(-lineage) %>%
  inner_join(select(cluster_meta, name, 'mean_pt'=pseudotime_ranks, lineage), by=c('RNA_snn_res.20'='name')) %>%
  arrange(mean_pt) %>%
  mutate(name=factor(RNA_snn_res.20, levels=unique(.$RNA_snn_res.20)), lineage=factor(lineage, levels=rev(names(lineage_colors))))


p_ex <- ggplot(ex_df, aes(name, gene, fill=expr)) +
  geom_tile() +
  facet_grid(gene_group~lineage, scales='free', space='free') +
  # scale_fill_gradientn(colors=pals::ocean.curl(100)) +
  scale_fill_gradientn(colors=greys(1)) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0.03, b = 0.03, r = 0.01, l = 0.01, unit='cm')
  )

p_mot <- ggplot(motif_df, aes(name, tf, fill=expr)) +
  geom_tile() +
  facet_grid(motif_group~lineage, scales='free', space='free') +
  # scale_fill_gradientn(colors=pals::ocean.amp(100))
  scale_fill_gradientn(colors=ylorbr(0.8)) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0.03, b = 0.03, r = 0.01, l = 0.01, unit='cm')
  )

p_pt <- ggplot(ex_df, aes(name, gene, fill=pseudotime_ranks)) +
  geom_tile() +
  facet_grid(~lineage, scales='free_x', space='free_x') +
  scale_fill_viridis(option='magma', direction=-1) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0.01, b = 0.01, r = 0.01, l = 0.01, unit='cm')
  )

p_line <- ggplot(meta, aes(x=name, fill=factor(line))) +
  geom_bar(position='fill') +
  facet_grid(~lineage, scales='free', space='free') +
  scale_fill_manual(values=line_colors) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0.01, b = 0.01, r = 0.01, l = 0.01, unit='cm')
  )

p_age <- ggplot(meta, aes(x=name, fill=factor(age))) +
  geom_bar(position='fill') +
  facet_grid(~lineage, scales='free', space='free') +
  scale_fill_manual(values=age_colors) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0.01, b = 0.01, r = 0.01, l = 0.01, unit='cm')
  )

p_lin <- ggplot(ex_df, aes(x=name, y=1, fill=factor(lineage))) +
  geom_tile() +
  facet_grid(~lineage, scales='free', space='free') +
  scale_fill_manual(values=lineage_colors) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0.01, b = 0.01, r = 0.01, l = 0.01, unit='cm')
  )

p_lin / p_pt / p_line / p_age / p_ex / p_mot + plot_layout(heights=c(0.8,0.8,1.5,1.5,20,20)) &
  no_legend() &
  scale_x_discrete(expand=c(0,0)) &
  scale_y_discrete(expand=c(0,0)) &
  theme(
    strip.text = element_blank(),
    panel.spacing.x = unit(0.02, 'cm'),
    panel.spacing.y = unit(0.02, 'cm'),
    panel.border = element_rect(color='grey', size=0.1, fill=NA)
  )
ggsave('plots/paper/fig1/expr_motif_heatmap.pdf', height=8, width=8, unit='cm')

p_ex / p_mot & theme(strip.text = element_blank(), axis.text.y=element_text()) & no_legend()
ggsave('plots/paper/fig1/expr_motif_heatmap_labels.pdf', height=40, width=30, unit='cm')

print_scale(ylorbr(0.8))
ggsave('plots/paper/fig1/ylorbr.pdf')
print_scale(greys(1))
ggsave('plots/paper/fig1/greys.pdf')



#### Tree construction ####
lineage_meta <- lineage_tree %>% as_tibble() %>% select(name, lineage)
# lineage_graph <- lineage_graph %>%
#   inner_join(select(hc_meta, cluster, UMAP_1, UMAP_2), by=c('name'='cluster')) %>%
#   inner_join(lineage_meta)

lineage_tree <- lineage_tree %N>%
  mutate(coords_x=case_when(
    lineage=='ge' ~ 0,
    lineage=='ctx' ~ -1,
    TRUE ~ coords_x
  ))

ggraph(lineage_graph, x=UMAP_1, y=UMAP_2) +
  geom_edge_link(alpha=0.2) +
  geom_node_point(aes(fill=to_ge_ranks), size=5, shape=21, stroke=0.2) +
  scale_fill_gradientn(colors=pals::ocean.deep(100)) +
  theme_void() + no_legend()
ggsave('plots/paper/fig1/graph_abstraction.png', width=5, height=3)


ggraph(lineage_graph, x=UMAP_1, y=UMAP_2) +
  geom_edge_link(alpha=0.2) +
  geom_node_point(size=2, shape=21, stroke=0.2, fill='gray') +
  scale_fill_gradientn(colors=pals::ocean.deep(100)) +
  theme_void() + no_legend()
ggsave('plots/paper/fig1/graph_basic.png', width=5, height=3)


ggraph(lineage_graph, x=UMAP_1, y=UMAP_2) +
  geom_edge_link(alpha=0.2) +
  geom_node_point(aes(fill=lineage), size=5, shape=21, stroke=0.2) +
  scale_fill_manual(values=lineage_colors) +
  theme_void() + no_legend()
ggsave('plots/paper/fig1/graph_lineages.png', width=5, height=3)


ggraph(lineage_graph, x=FR1, y=FR2) +
  geom_edge_link(alpha=0.2, width=0.5) +
  geom_node_point(aes(fill=lineage), size=3, shape=21, stroke=0.2) +
  scale_fill_manual(values=lineage_colors) +
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_void() + no_legend()
ggsave('plots/paper/fig1/tree_lineages.png', width=4.3, height=3.7)

pt_size <- 4
p_ctx <- ggraph(lineage_tree, x=pseudotime_align, y=coords_x) +
  geom_edge_link(alpha=0.2) +
  geom_node_point(aes(fill=to_ctx_ranks), size=pt_size, shape=21, stroke=0.2) +
  scale_fill_gradientn(colors=pals::ocean.deep(100)) +
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_void() + no_legend()

p_ge <- ggraph(lineage_tree, x=pseudotime_align, y=coords_x) +
  geom_edge_link(alpha=0.2) +
  geom_node_point(aes(fill=to_ge_ranks), size=pt_size, shape=21, stroke=0.2) +
  scale_fill_gradientn(colors=pals::ocean.deep(100)) +
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_void() + no_legend() + theme(plot.margin = unit(c(0,1,0,0), 'cm'))

p_nt <- ggraph(lineage_tree, x=pseudotime_align, y=coords_x) +
  geom_edge_link(alpha=0.2) +
  geom_node_point(aes(fill=to_nt_ranks), size=pt_size, shape=21, stroke=0.2) +
  scale_fill_gradientn(colors=pals::ocean.deep(100)) +
  scale_y_continuous(expand=c(0.1,0.1)) +
  theme_void() + no_legend() + theme(plot.margin = unit(c(0,1,0,0), 'cm'))

p_nt + p_ge + p_ctx
ggsave('plots/paper/fig1/tree_transition.png', width=7.5, height=1.5)

print_scale(pals::ocean.deep(100))
ggsave('plots/paper/fig1/ocean.pdf')


#### Select features to plot ####
peaks_near_gene <- find_peaks_near_genes(peak_ranges, gene_annot, distance = 50000, only_tss = FALSE)
peaks2gene <- aggregate.Matrix(t(peaks_near_gene), groupings=colnames(peaks_near_gene))
peaks2gene_df <- peaks2gene %>% summary()
peaks2gene_df$gene <- rownames(peaks2gene)[peaks2gene_df$i]
peaks2gene_df$feature <- colnames(peaks2gene)[peaks2gene_df$j]
peaks2gene_df <- peaks2gene_df %>%
  as_tibble() %>%
  select(gene, feature)


telen_peaks <- de_lin_list$peaks %>%
  rename('feature'='gene') %>%
  mutate(expr_ratio=prcex_self/prcex_other) %>%
  filter(group=='telencephalon', prcex_self>10) %>%
  inner_join(peaks2gene_df) %>%
  top_n(80, expr_ratio)






#### Plot on graph ####
cluster_summaries <- read_rds('data/RNA_ATAC/RNA_ATAC_res20_summaries.rds')


##### RNA ####
gene_expr_df <- cluster_summaries$RNA[, pick_genes] %>%
  {colnames(.) <- str_replace(colnames(.), '-', '_');.} %>%
  as_tibble(rownames='name')

perc_expr_df <- cluster_summaries$perc_expr[, pick_genes] %>%
  {colnames(.) <- str_replace(colnames(.), '-', '_');.} %>%
  {colnames(.) <- paste0(colnames(.), '_perc');.} %>%
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>%
  left_join(gene_expr_df) %>%
  left_join(perc_expr_df)

plots <- map(str_replace(pick_genes, '-', '_'), function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes_string(fill=gene, size=paste0(gene, '_perc')), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.5, 3)) +
    scale_fill_gradientn(colors=gyylgnbu()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend() + ggtitle(gene)
})
p_ex <- wrap_plots(plots)
p_ex

##### gene acivity ####
pick_genes_motifs <- motif2tf %>% filter(tf%in%pick_genes) %>% pull(motif)
gene_motifs_df <- cluster_summaries$chromvar[, pick_genes_motifs] %>%
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>%
  left_join(gene_motifs_df)

plots <- map(str_replace(pick_genes_motifs, '-', '_'), function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes_string(fill=gene), size=2.5, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=gyylorbr()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend() + ggtitle(filter(motif2tf, motif==gene)$name[1])
})
p_act <- wrap_plots(plots)




##### Peaks ####
pick_peaks <- unique(telen_peaks$feature)
peak_acc_df <- cluster_summaries$peaks[, pick_peaks] %>%
  {colnames(.) <- str_replace_all(colnames(.), '-', '_');.} %>%
  as_tibble(rownames='name')

perc_acc_df <- cluster_summaries$perc_acc[, pick_peaks] %>%
  {colnames(.) <- str_replace_all(colnames(.), '-', '_');.} %>%
  {colnames(.) <- paste0(colnames(.), '_perc');.} %>%
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>%
  left_join(peak_acc_df) %>%
  left_join(perc_acc_df)

plots <- map(str_replace_all(pick_peaks, '-', '_'), function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    # geom_edge_link(alpha=0.1) +
    geom_node_point(aes_string(fill=gene, size=paste0(gene, '_perc')), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.5, 3)) +
    scale_fill_gradientn(colors=gyylgnbu()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend() + ggtitle(gene)
})
p_acc <- wrap_plots(plots)
p_acc


#### Final feature plots on graph ####
cluster_summaries <- read_rds('data/RNA_ATAC/RNA_ATAC_res20_summaries.rds')

pick_genes <- c('EGR1', 'FOXG1', 'NFIA', 'NFIB', 'DLX5', 'ZIC1', 'ZIC2', 'LHX5')

gene_expr_df <- cluster_summaries$RNA[, pick_genes] %>%
  {colnames(.) <- str_replace(colnames(.), '-', '_');.} %>%
  as_tibble(rownames='name')

perc_expr_df <- cluster_summaries$perc_expr[, pick_genes] %>%
  {colnames(.) <- str_replace(colnames(.), '-', '_');.} %>%
  {colnames(.) <- paste0(colnames(.), '_perc');.} %>%
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>%
  left_join(gene_expr_df) %>%
  left_join(perc_expr_df)

plots <- map(str_replace(pick_genes, '-', '_'), function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes_string(fill=gene, size=paste0(gene, '_perc')), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.5, 3)) +
    scale_fill_gradientn(colors=gyylgnbu()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend() + ggtitle(gene)
})
p_ex <- wrap_plots(plots)
p_ex


pick_genes_motifs <- motif2tf %>% filter(tf%in%pick_genes) %>% pull(motif)
gene_motifs_df <- cluster_summaries$chromvar[, pick_genes_motifs] %>%
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>%
  left_join(gene_motifs_df)

plots <- map(str_replace(pick_genes_motifs, '-', '_'), function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes_string(fill=gene), size=2.5, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=gyylorbr()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend() + ggtitle(filter(motif2tf, motif==gene)$name[1])
})
p_act <- wrap_plots(plots)

p_ex | p_act
# ggsave('plots/paper/fig1/pick_expr_motifs_graph.png', width=16, height=8)


#### Pt bins ####
var_feats <- FindVariableFeatures(rnatac, nfeatures = 1000) %>% VariableFeatures()


rnatac$pt_bin <- cut(rnatac$pseudotime_ranks, breaks = 20, labels = 1:20)
rnatac$stage_manual <- case_when(
  rnatac$pseudotime_ranks < 0.1 ~ 'iPSC',
  rnatac$pseudotime_ranks < 0.4 ~ 'nect_nepi',
  rnatac$pseudotime_ranks < 0.9 ~ 'npc',
  TRUE ~ 'neuron'
)

pt_summary <- GetAssayData(rnatac, assay='RNA', slot='data')[var_feats, ] %>% t() %>% scale() %>%
  aggregate_matrix(groups=rnatac$pt_bin)

pt_clust <- pt_summary %>% t() %>% sparse_cov() %>% {1-.$cor} %>% as.dist() %>% hclust()
pt_clust$order <- 1:20

plot(pt_clust)


pt_clusters <- pt_clust %>% cutree(h=1)
rnatac$pt_clusters <- pt_clusters[rnatac$pt_bin]

rnatac$stage_manual_hclust <- case_when(
  as.numeric(as.character(rnatac$pt_bin)) <= 2 ~ 'iPSC',
  as.numeric(as.character(rnatac$pt_bin)) <= 7 ~ 'nect',
  as.numeric(as.character(rnatac$pt_bin)) <= 11 ~ 'nepi',
  as.numeric(as.character(rnatac$pt_bin)) <= 17 ~ 'npc',
  TRUE ~ 'neuron'
)

atac_meta <- rnatac@meta.data[ps_order, c('stage_manual_hclust'), drop=F]
atac_meta <- rownames_to_column(atac_meta, 'pseudocell')
rownames(atac_meta) <- colnames(atac)
atac <- AddMetaData(atac, atac_meta)

dim_plot(atac, group.by=c('stage', 'stage_manual_hclust', 'lineage'))
dim_plot(rnatac, group.by=c('stage', 'stage_manual_hclust', 'lineage'))


p1 <- dim_plot(rnatac, group.by=c('pt_bin', 'stage_manual', 'stage_manual_hclust'))
p2 <- feature_plot(rnatac, features=c('POU5F1', 'KRT8', 'HES1', 'NOTCH1', 'HES5', 'VIM', 'SOX3', 'DCX'), order=T)
p1 / p2 + plot_layout(heights=c(1,2))



stage_peak_de <- de(rnatac, assay = 'peaks', groups = 'stage_manual_hclust')
stage_peak_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_acc_stage_coarse_hclust.tsv')



#### DE between stages (again) ####
expr_stage_de <- de(rnatac, 'stage_manual_hclust', assay='RNA')
motif_stage_de <- de(rnatac, 'stage_manual_hclust', assay='chromvar')
act_stage_de <- de(rnatac, 'stage_manual_hclust', assay='gene_activity')

expr_stage_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_expr_stage_hclust.tsv')
motif_stage_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_motif_stage_hclust.tsv')
act_stage_de %>% write_tsv('data/RNA_ATAC/diff_expression/diff_act_stage_hclust.tsv')

expr_stage_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_expr_stage_hclust.tsv')
motif_stage_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_motif_stage_hclust.tsv')
act_stage_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_act_stage_hclust.tsv')


expr_stage_de


#### Get stage specific peaks ####

stage_peak_de <- read_tsv('data/RNA_ATAC/diff_expression/diff_acc_stage_hclust.tsv')

stage_peak_sig <- stage_peak_de %>%
  mutate(
    expr_frac=(prcex_self/prcex_other),
    chrom=str_replace(feature, '(chr[\\d\\w]+)-\\d+-\\d+', '\\1')
  ) %>%
  group_by(feature) %>%
  filter(expr_frac==max(expr_frac)) %>%
  group_by(group) %>%
  filter(prcex_self>10, expr_frac>1.5) %>%
  top_n(5000, expr_frac)


ggplot(stage_peak_sig, aes(chrom)) +
  geom_bar() +
  facet_wrap(~group)

vulcano_plot(stage_peak_sig)


stage_de <- bind_rows('RNA'=expr_stage_de, 'chromvar'=motif_stage_de, 'gene_activity'=act_stage_de, .id='assay')
stage_de <- stage_de %>% filter(padj<0.01)

stage_peak_sig$assay <- 'peaks'
stage_de_all <- bind_rows(stage_peak_sig, stage_de)

stage_de_all %>% write_tsv('data/RNA_ATAC/diff_expression/all_stage_hclust_de.tsv')


### Cluster peaks (for each group separately) ####
stage_split <- stage_peak_sig %>% group_by(group) %>% group_split()
names(stage_split) <- unique(stage_peak_sig$group)
stage_summary <- map_par(stage_split, function(df){
  group_mask <- rnatac$stage_manual_hclust%in%df$group
  peaks <- df %>% pull(feature) %>% unique()
  peak_summary <- GetAssayData(rnatac, assay='peaks', slot='data')[peaks, group_mask]
  peak_clust <- peak_summary %>% t() %>% sparse_cov() %>% {1-.$cor} %>% as.dist() %>% hclust()
  return(list(
    agg = peak_summary,
    clust = peak_clust
  ))
}, parallel=TRUE)

map(stage_summary, function(x) {cutree(x$clust, h=1) %>% max()})



#### Accessible genome size for each pt bin ####

peak_mat <- GetAssayData(rnatac, assay='peaks')
accessible_genome_size <- map_par(1:ncol(peak_mat), function(i){
  cell_peaks <- peak_mat[,i]
  cell_peaks <- names(cell_peaks[cell_peaks>0])
  peak_width <- StringToGRanges(cell_peaks) %>% reduce() %>% width() %>% sum()
  return(peak_width)
}, parallel=T)

acc_genome <- read_rds('data/RNA_ATAC/acc_genome_size.rds')
acc_genome <- acc_genome %>% enframe('cell', 'acc_bp')


#### Plot heatmap ####
stage_cluster_df <- stage_summary %>% map_dfr(function(x) enframe(cutree(x$clust, h=1)[x$clust$order]), .id='group') %>%
  mutate(group_clust = paste0(group, '_', value)) %>%
  group_by(group) %>%
  mutate(order=as.numeric(factor(group_clust, levels=unique(.$group_clust))))

pt_summary <- GetAssayData(rnatac, assay='peaks', slot='data')[stage_cluster_df$name, ] %>% t() %>%
  aggregate_matrix(groups=as.character(rnatac$pt_bin))

peak_clust_summary <- scale(pt_summary, center=F) %>% t() %>% aggregate_matrix(groups=stage_cluster_df$group_clust)

cluster_order <- peak_clust_summary %>% t() %>% sparse_cov() %>% {1-.$cor} %>% as.dist() %>% hclust(method='ward.D2') %>% {.$label[.$order]}



clust_meta <- stage_cluster_df %>%
  distinct(group, group_clust, order)
peak_clust_summary_df <- peak_clust_summary %>%
  as_tibble(rownames='group_clust') %>%
  pivot_longer(!group_clust) %>%
  inner_join(clust_meta) %>%
  mutate(
    name=factor(name, levels=sort(unique(as.numeric(name)))),
    group=factor(group, levels=c('iPSC', 'nect', 'nepi', 'npc', 'neuron'))
  ) %>%
  mutate(group_clust=factor(group_clust, levels=cluster_order)) %>%
  arrange(group, group_clust) %>%
  mutate(group_clust=factor(group_clust, levels=rev(unique(.$group_clust)))) %>%
  group_by(group_clust) %>%
  mutate(value=scale01(value))

pt_meta <- rnatac@meta.data %>%
  as_tibble(rownames='cell') %>%
  inner_join(acc_genome)

p1 <- ggplot(peak_clust_summary_df, aes(name, group_clust, fill=value)) +
  geom_tile() +
  facet_grid(group~., scales='free', space='free') +
  scale_fill_gradientn(colors=greys(1.5)) +
  article_text() +
  no_margin() +
  no_y_text() + no_legend() +
  no_x_text() +
  theme(
   panel.border = element_blank(),
   strip.text.y = element_blank(),
   panel.spacing = unit(0.05, 'lines')
  ) +
  labs(x='pseudotime bins', y='peak clusters', fill='accessibility')

p2 <- ggplot(pt_meta, aes(pt_bin, 1, fill=stage_manual_hclust)) +
  geom_tile() +
  scale_fill_manual(values=stage_colors) +
  theme_void() + no_legend()

p3 <- ggplot(pt_meta, aes(pt_bin, fill=factor(age))) +
  geom_bar(position='fill') +
  scale_fill_manual(values=age_colors) +
  theme_void() + no_legend() + no_margin()

p4 <- ggplot(pt_meta, aes(pt_bin, acc_bp/1000000)) +
  geom_violin(fill='darkgrey', color='black', size=0.1) +
  geom_boxplot(fill='white', color='black', width=0.1, outlier.shape=NA, size=0.1) +
  theme_rangeframe() + scale_axis_rangeframe() +
  scale_y_continuous(breaks=c(0,100)) +
  article_text() +
  no_legend() + no_margin() +
  no_x_text() +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y='Total accessible\nregion [Mb]')

p4 <- ggplot(pt_meta, aes(pt_bin, acc_bp/1000000)) +
  geom_violin(fill='darkgrey', color='black', size=0.1) +
  geom_boxplot(fill='white', color='black', width=0.1, outlier.shape=NA, size=0.1) +
  theme_rangeframe() + scale_axis_rangeframe() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10), breaks=c(10,50,100)) +
  article_text() +
  no_legend() + no_margin() +
  no_x_text() +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y='Total accessible\nregion [Mb]')

p4

p4 / p2 / p3 / p1 + plot_layout(heights=c(8,1,2,30))
ggsave('plots/paper/fig1/et_acc_chrom_heatmap.pdf', height=5, width=4.5, units='cm')


pdf('plots/paper/fig1/peak_cluster_dend.pdf')
plot(pt_clust, hang = -1, cex = 0.6)
dev.off()






#### Write peaks for GREAT ####
stage_cluster_peaks_list <- stage_cluster_df %>%
  group_by(group_clust) %>% group_split()


background_peaks <- StringToGRanges(rownames(rnatac@assays$peaks))
names(background_peaks) <- rownames(rnatac@assays$peaks)
background_peaks$name <- rownames(rnatac@assays$peaks)
export.bed(background_peaks, 'data/RNA_ATAC/GREAT/background_peaks.bed')

stage_cluster_peaks_list <- stage_cluster_df %>%
  group_by(group_clust) %>% group_split()
map(stage_cluster_peaks_list, function(df){
  peaks <- StringToGRanges(df$name)
  names(peaks) <- df$name
  peaks$name <- df$name
  filename <- paste0('data/RNA_ATAC/GREAT/', unique(df$group_clust), '_peaks.bed')
  export.bed(peaks, filename)
})


stage_peaks_list <- stage_cluster_df %>%
  group_by(group) %>% group_split()
map(stage_peaks_list, function(df){
  peaks <- StringToGRanges(df$name)
  names(peaks) <- df$name
  peaks$name <- df$name
  filename <- paste0('data/RNA_ATAC/GREAT/', unique(df$group), '_peaks.bed')
  export.bed(peaks, filename)
})



### Read GREAT results

great_results <- list.files('data/RNA_ATAC/GREAT/results/', pattern='.+_all.tsv', full.names = T) %>% {names(.)<-.;.} %>% map_dfr(function(x) read_tsv(x, skip=3), .id='file') %>%
  rename('Ontology'=2) %>%
  mutate(stage=str_replace(file, '.+results//(\\w+)_all.tsv', '\\1'))

great_results %>% select(-file) %>% filter(HyperFdrQ<0.05) %>% write_tsv('data/RNA_ATAC/GREAT/results/all_stages.tsv')

great_results_sig <- great_results %>%
  filter(HyperFdrQ<1e-2, ObsRegions>30)


great_results_top <- great_results_sig %>%
  group_by(stage) %>%
  top_n(20, RegionFoldEnrich) %>%
  arrange(RegionFoldEnrich) %>%
  mutate(Desc=factor(Desc, levels=unique(.$Desc)))

ggplot(great_results_top, aes(RegionFoldEnrich, Desc, label=Desc)) +
  geom_bar(stat='identity', fill='#ebedef', color='#2e4053', size=0.1) +
  geom_text(aes(x=0), hjust=0, color='black', nudge_x=0.1, size=10/ggplot2::.pt) +
  scale_x_continuous(expand=c(0,0)) +
  facet_wrap(~stage, scales='free') +
  theme_rangeframe() + scale_axis_rangeframe() +
  no_y_text() + theme(axis.line.y = element_blank()) +
  labs(y='', x='Region fold enrichment')



great_results_top_list <- great_results_sig %>%
  group_by(stage) %>%
  top_n(17, RegionFoldEnrich) %>%
  group_split() %>%
  map(function(df){
    df %>%
      arrange(RegionFoldEnrich) %>%
      mutate(Desc=factor(Desc, levels=unique(.$Desc))) %>%
      return()
  })

plots <- map(great_results_top_list, function(df){
  ggplot(df, aes(RegionFoldEnrich, Desc, label=Desc)) +
    geom_bar(stat='identity', fill='#ebedef', color='#2e4053', size=0.1) +
    geom_text(aes(x=0), hjust=0, color='black', nudge_x=0.1, size=5/ggplot2::.pt) +
    scale_x_continuous(expand=c(0,0)) +
    facet_wrap(~stage, scales='free') +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_y_text() + theme(axis.line.y = element_blank()) +
    no_margin() +
    labs(y='', x='Region fold enrichment')
})

wrap_plots(plots, ncol=2)
ggsave('plots/paper/supp_et/atac_great_enrich_bar.pdf', width=9.7, height=13, units='cm')




#### Tracks ####
atac_cov <- subset(atac, lineage%in%c('nt', 'telencephalon', 'early', 'ge', 'ctx'))
atac_cov$lineage <- factor(atac_cov$lineage, levels=c('early', 'nt', 'telencephalon', 'ge', 'ctx'))


atac_meta <- rnatac@meta.data[ps_order, c('RNA_snn_res.20', 'lineage', 'lineage_coarse', 'stage_manual_hclust')]
atac_meta <- rownames_to_column(atac_meta, 'pseudocell')
rownames(atac_meta) <- colnames(atac)
atac_cov <- AddMetaData(atac_cov, atac_meta)
atac_cov$stage_manual_hclust <- factor(atac_cov$stage_manual_hclust, levels=c('iPSC', 'nect', 'nepi', 'npc', 'neuron'))


#### POU5F1 ####
p_cov <- CoveragePlot(
  atac_cov,
  group.by = 'stage_manual_hclust',
  region = 'POU5F1',
  extend.upstream = 4000,
  extend.downstream = 6000,
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=stage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

region <- FindRegion(
  object = atac_cov,
  region = 'POU5F1',
  assay = 'peaks',
  extend.upstream = 4000,
  extend.downstream = 6000,
)

p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void()

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()


p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))
ggsave('plots/paper/fig1/pou5f1_tracks.pdf', unit='cm')


#### FGF8 ####
goi <- 'FGF8'

annot <- gene_annot[gene_annot$gene_name==goi]
strand <- mode(as.character(strand(annot)))
region <- GRanges(
  seqnames = as.character(seqnames(annot))[[1]],
  ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
  strand = strand
)
region <- resize(region, width=0, fix='end')
region <- Extend(region, upstream = 500, downstream=12000)


p_cov <- CoveragePlot(
  atac_cov,
  group.by = 'stage_manual_hclust',
  region = region,
  extend.upstream = 0,
  extend.downstream = 0,
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=stage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

region <- FindRegion(
  object = atac_cov,
  region = region,
  assay = 'peaks',
  extend.upstream = 0,
  extend.downstream = 0,
)

p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void()

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()


p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))

ggsave('plots/paper/fig1/fgf8_tracks.pdf', unit='cm')


#### VIM ####
goi <- 'VIM'

annot <- gene_annot[gene_annot$gene_name==goi]
strand <- mode(as.character(strand(annot)))
region <- GRanges(
  seqnames = as.character(seqnames(annot))[[1]],
  ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
  strand = strand
)
region <- resize(region, width=0, fix='end')
region <- Extend(region, upstream = 6500, downstream=3000)


p_cov <- CoveragePlot(
  atac_cov,
  group.by = 'stage_manual_hclust',
  region = region,
  extend.upstream = 0,
  extend.downstream = 0,
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=stage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

region <- FindRegion(
  object = atac_cov,
  region = region,
  assay = 'peaks',
  extend.upstream = 0,
  extend.downstream = 0,
)

p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void()

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()


p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))


ggsave('plots/paper/fig1/vim_tracks.pdf', unit='cm')


#### STMN4 ####
goi <- 'STMN4'

annot <- gene_annot[gene_annot$gene_name==goi]
strand <- mode(as.character(strand(annot)))
region <- GRanges(
  seqnames = as.character(seqnames(annot))[[1]],
  ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
  strand = strand
)
region <- resize(region, width=0, fix='start')
region <- Extend(region, upstream = 8000, downstream=5000)


p_cov <- CoveragePlot(
  atac_cov,
  group.by = 'stage_manual_hclust',
  region = region,
  extend.upstream = 0,
  extend.downstream = 0,
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=stage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

region <- FindRegion(
  object = atac_cov,
  region = region,
  assay = 'peaks',
  extend.upstream = 0,
  extend.downstream = 0,
)

p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void()

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()


p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))


ggsave('plots/paper/fig1/stmn4_tracks.pdf', unit='cm')



#### DCX ####
goi <- 'DCX'

annot <- gene_annot[gene_annot$gene_name==goi]
strand <- mode(as.character(strand(annot)))
region <- GRanges(
  seqnames = as.character(seqnames(annot))[[1]],
  ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
  strand = strand
)
region <- resize(region, width=0, fix='start')
region <- Extend(region, upstream = -1200, downstream=6000)


p_cov <- CoveragePlot(
  atac_cov,
  group.by = 'stage_manual_hclust',
  region = region,
  extend.upstream = 0,
  extend.downstream = 0,
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=stage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

region <- FindRegion(
  object = atac_cov,
  region = region,
  assay = 'peaks',
  extend.upstream = 0,
  extend.downstream = 0,
)

p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void()

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()


p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))


ggsave('plots/paper/fig1/dcx_tracks.pdf', unit='cm')



#### BCL11A ####
goi <- 'EMX1'

annot <- gene_annot[gene_annot$gene_name==goi]
strand <- mode(as.character(strand(annot)))
region <- GRanges(
  seqnames = as.character(seqnames(annot))[[1]],
  ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
  strand = strand
)
region <- resize(region, width=0, fix='start')
region <- Extend(region, upstream = 20000, downstream=5000)


p_cov <- CoveragePlot(
  atac_cov,
  group.by = 'lineage',
  region = region,
  extend.upstream = 0,
  extend.downstream = 0,
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=stage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

region <- FindRegion(
  object = atac_cov,
  region = region,
  assay = 'peaks',
  extend.upstream = 0,
  extend.downstream = 0,
)

p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void()

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()


p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))


ggsave('plots/paper/fig1/dcx_tracks.pdf', unit='cm')




# [1] "ATP6AP2,BMP4,EN1,FGF8,FOXA1,GBX2,GLI2,GLI3,GSC,HES1,HES3,PAX6,PAX7,PSEN1,PTCH1,SHH,SOX17,SSBP3,TCTN1"



#### FOXG1 ####
top_act <- act_de %>%
  filter(group=='organoid') %>%
  top_n(100, abs(fc))


p_cov <- CoveragePlot(
  atac_cov,
  group.by = 'stage',
  region = 'FOXG1',
  extend.upstream = 10000,
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=stage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

region <- FindRegion(
  object = atac_cov,
  region = 'FOXG1',
  assay = 'peaks',
  extend.upstream = 10000
)


p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void()

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()


p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))
ggsave('plots/paper/fig1/foxg1_tracks.pdf', unit='cm')

PeakPlot(
  atac_cov,
  region = region
)
ggsave('plots/paper/fig1/foxg1_peaks.pdf', unit='cm', height=3)



#### KRT8 ####
top_act <- act_de %>%
  filter(group=='organoid') %>%
  top_n(30, -fc)

CoveragePlot(atac_cov, group.by='stage', region = 'KRT8', extend.upstream = 10000)


p_cov <- CoveragePlot(
  atac_cov,
  group.by = 'stage',
  region = 'KRT8',
  extend.upstream = 5000,
  annotation = FALSE,
  peaks = FALSE
) & scale_fill_manual(values=stage_colors) &
  theme_void() & no_legend() & theme(strip.text = element_blank())

region <- FindRegion(
  object = atac_cov,
  region = 'KRT8',
  assay = 'peaks',
  extend.upstream = 5000
)


p_an <- AnnotationPlot(
  atac_cov,
  region = region
) + theme_void()

p_pe <- PeakPlot(
  atac_cov,
  region = region
) + theme_void() + no_legend()

p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))
ggsave('plots/paper/fig1/krt8_tracks.pdf', unit='cm')

PeakPlot(
  atac_cov,
  region = region
)
ggsave('plots/paper/fig1/krt8_peaks.pdf', unit='cm', height=3)




#### DLX5 ####
CoveragePlot(atac_cov, group.by='lineage', region = 'DLX5', extend.upstream = 10000)

region <- FindRegion(
  object = atac_cov,
  region = 'DLX5',
  assay = 'peaks',
  extend.upstream = 9000
)

region <- resize(region, width = 9500, fix='start')

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

p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))
ggsave('plots/paper/fig1/dlx5_tracks.pdf', unit='cm')

PeakPlot(
  atac_cov,
  region = region
)
ggsave('plots/paper/fig1/dlx5_peaks.pdf', unit='cm', height=3)





#### NFIA ####
CoveragePlot(atac_cov, group.by='lineage', region = 'NFIA', extend.upstream = 10000)

region <- FindRegion(
  object = atac_cov,
  region = 'NFIA',
  assay = 'peaks',
  extend.upstream = 9000
)

region <- resize(region, width = 10000, fix='start')

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

p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))
ggsave('plots/paper/fig1/nfia_tracks.pdf', unit='cm')

PeakPlot(
  atac_cov,
  region = region
)
ggsave('plots/paper/fig1/nfia_peaks.pdf', unit='cm', height=3)




#### FOXG1 ####
CoveragePlot(atac_cov, group.by='lineage', region = 'FOXG1', extend.upstream = 10000)

region <- FindRegion(
  object = atac_cov,
  region = 'FOXG1',
  assay = 'peaks',
  extend.upstream = 5500
)

region <- resize(region, width = 6000, fix='start')

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

p_cov / p_an / p_pe + plot_layout(heights=c(10,1,1)) & theme(text = element_text(size=5))
ggsave('plots/paper/fig1/foxg1_tracks.pdf', unit='cm')

PeakPlot(
  atac_cov,
  region = region
)
ggsave('plots/paper/fig1/foxg1_peaks.pdf', unit='cm', height=3)























