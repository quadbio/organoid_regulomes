source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')

setwd('~/projects/early/')

library(SeuratDisk)


#### Read Stuff ####
rnatac <- read_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2.1_srt.rds')
rna <- read_rds('data/RNA/RNA_all_merged/RNA_all_merged_v2.1_srt.rds')
atac <- read_rds('data/ATAC/ATAC_all_merged/ATAC_all_merged_v3.1_srt.rds')

rnatac_res20_summaries <- read_rds('data/RNA_ATAC/RNA_ATAC_res20_summaries.rds')
rnatac_res20_meta <- read_tsv('data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv') %>% 
  mutate(cluster=as.character(cluster))

motif2tf <- read_tsv('data/tf_matches/motif2tf_early_all.tsv') 

lineage_graph <- read_rds('data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')
lineage_meta <- read_tsv('data/RNA_ATAC/lineages/cellrank_graph_res20_meta.tsv')


lineage_meta <- as_tibble(lineage_graph) %>% 
  pivot_longer(cols=c(to_ctx, to_ge, to_nt), names_to='to_fate')
ggplot(lineage_meta, aes(pseudotime_ranks, value, fill=lineage)) +
  geom_point(size=0.7, shape=21, stroke=0.1, color='white') +
  facet_grid(~to_fate) +
  scale_fill_manual(values=lineage_colors) +
  theme_rangeframe() + scale_axis_rangeframe() +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  article_text() +
  no_legend() + 
  theme(
    strip.text = element_blank()
  ) +
  labs(x = 'Velocity pseudotime', y='Transition probability')
ggsave('plots/paper/supp_et/lineage_probs_raw_pt_scatter.pdf', width=6.8, height=3.4, units='cm')


meta_df <- lineage_meta %>% column_to_rownames('name')

cluster_srt <- CreateSeuratObject(t(rnatac_res20_summaries$RNA[rownames(meta_df), ]))
cluster_srt <- AddMetaData(cluster_srt, meta_df)
cluster_srt[['fr']] <- CreateDimReducObject(as.matrix(meta_df[colnames(cluster_srt), c('FR1', 'FR2')]), key='FR_')
cluster_srt[['umap']] <- CreateDimReducObject(as.matrix(meta_df[colnames(cluster_srt), c('UMAP_1', 'UMAP_2')]), key='UMAP_')
cluster_srt[['chromvar']] <- CreateAssayObject(t(rnatac_res20_summaries$chromvar[rownames(meta_df), ]))
cluster_srt[['gene_activity']] <- CreateAssayObject(t(rnatac_res20_summaries$gene_activity[rownames(meta_df), ]))
cluster_srt[['peaks']] <- CreateAssayObject(t(rnatac_res20_summaries$peaks[rownames(meta_df), ]))
cluster_srt[['module_score']] <- CreateAssayObject(t(rnatac_res20_summaries$module_score[rownames(meta_df), ]))
cluster_srt[['module_neg_score']] <- CreateAssayObject(t(rnatac_res20_summaries$module_neg_score[rownames(meta_df), ]))
cluster_srt[['module_peaks_chromvar']] <- CreateAssayObject(t(rnatac_res20_summaries$module_peaks_chromvar[rownames(meta_df), ]))
cluster_srt[['module_neg_peaks_chromvar']] <- CreateAssayObject(t(rnatac_res20_summaries$module_neg_peaks_chromvar[rownames(meta_df), ]))
cluster_srt[['module_all_peaks_chromvar']] <- CreateAssayObject(t(rnatac_res20_summaries$module_all_peaks_chromvar[rownames(meta_df), ]))
cluster_srt[['perc_expr']] <- CreateAssayObject(t(rnatac_res20_summaries$perc_expr[rownames(meta_df), ]))
cluster_srt[['perc_acc']] <- CreateAssayObject(t(rnatac_res20_summaries$perc_acc[rownames(meta_df), ]))

cluster_srt %>% write_rds('data/RNA_ATAC/RNA_ATAC_res20_srt.rds')

dim_plot(cluster_srt, reduction='fr', group.by='lineage')
feature_plot(cluster_srt, reduction='fr', features=c('NEUROD6', 'PAX6', 'SIX6', 'EOMES', 'HOPX', 'VIM', ''))


#### QC plots ####
rna_meta <- rna@meta.data %>% as_tibble(rownames='cell') %>% 
  arrange(age) %>% mutate(batch=factor(orig.ident, levels=unique(.$orig.ident)))

ggplot(rna_meta, aes(nFeature_RNA, fill=orig.ident)) +
  # geom_histogram() +
  geom_histogram(fill='darkgrey', color='black', size=0.2) +
  facet_grid(batch~., scales='free_y') +
  scale_x_continuous(limits=c(0, 9000), breaks=c(0, 4000, 8000)) +
  scale_y_continuous(breaks=scales::breaks_extended(3)) +
  # scale_fill_manual(values=age_colors) +
  theme_rangeframe() + scale_axis_rangeframe() +
  article_text() + 
  theme(
    strip.text.y = element_blank()
  ) +
  labs(y='Count', x='# detected features')
ggsave('plots/paper/supp_et/rna_qc_feature_hist.pdf', width=3.7, height=7, units='cm')



atac_meta <- atac@meta.data %>% as_tibble(rownames='cell') %>% 
  arrange(age) %>% mutate(batch=factor(batch, levels=unique(.$batch))) 
ggplot(atac_meta, aes(log10(peak_region_fragments), fill=factor(batch))) +
  # geom_histogram() +
  geom_histogram(fill='darkgrey', color='black', size=0.2) +
  facet_grid(batch~., scales='free_y') +
  # scale_x_continuous(limits=c(0, 9000), breaks=c(0, 4000, 8000)) +
  scale_y_continuous(breaks=scales::breaks_extended(3)) +
  # scale_fill_manual(values=age_colors) +
  theme_rangeframe() + scale_axis_rangeframe() +
  article_text() + 
  theme(
    strip.text.y = element_blank()
  ) +
  labs(y='Count', x=expression(log[10] ~'peak region fragments'))
ggsave('plots/paper/supp_et/atac_qc_peak_fracs_hist.pdf', width=3.7, height=7, units='cm')

# 
# TSSEnrichment(atac, assay = 'peaks')
# TSSPlot(atac, assay='peaks')


ggplot(atac_meta, aes(pct_reads_in_peaks, fill=factor(batch))) +
  # geom_histogram() +
  geom_histogram(fill='darkgrey', color='black', size=0.2) +
  facet_grid(batch~., scales='free_y') +
  # scale_x_continuous(limits=c(0, 9000), breaks=c(0, 4000, 8000)) +
  scale_y_continuous(breaks=scales::breaks_extended(3)) +
  # scale_fill_manual(values=age_colors) +
  theme_rangeframe() + scale_axis_rangeframe() +
  article_text() + 
  theme(
    strip.text.y = element_blank()
  ) +
  labs(y='Count', x='Percent reads in peaks')
ggsave('plots/paper/supp_et/atac_qc_reads_in_peaks_hist.pdf', width=3.7, height=7, units='cm')


#### Demux data ####
all_demux <- read_tsv('data/demux/ALL_demux_v2.tsv')

demux_rna_meta <- rna@meta.data %>% 
  as_tibble(rownames='cell') %>% 
  dplyr::select(-line, -p_singlet) %>% 
  left_join(all_demux) %>% 
  arrange(age) %>% mutate(orig.ident=factor(orig.ident, levels=unique(.$orig.ident)))

ggplot(demux_rna_meta) +
  geom_histogram(mapping=aes(llk1), fill=blue, alpha=0.5) +
  geom_histogram(mapping=aes(llk2), fill=redpink, alpha=0.5) +
  scale_y_continuous(breaks=scales::breaks_extended(3)) +
  theme_article() +
  # facet_grid(orig.ident~., scales='free') +
  theme_rangeframe() + scale_axis_rangeframe() +
  lims(x=c(-1000,0)) +
  article_text() +
  # theme(
  #   strip.text = element_blank()
  # ) + 
  labs(x='Log likelihood', y='Count')
  
ggsave('plots/paper/supp_et/rna_qc_demux_llk_hist.pdf', width=3.5, height=1.7, units='cm')


demux_atac_meta <- atac@meta.data %>% 
  as_tibble(rownames='cell') %>% 
  dplyr::select(-line, -p_singlet) %>% 
  left_join(all_demux) %>% 
  arrange(age) %>% mutate(batch=factor(batch, levels=unique(.$batch)))

ggplot(filter(demux_atac_meta, age!=4)) +
  geom_histogram(mapping=aes(llk1), fill=blue, alpha=0.5) +
  geom_histogram(mapping=aes(llk2), fill=redpink, alpha=0.5) +
  scale_y_continuous(breaks=scales::breaks_extended(3)) +
  theme_article() +
  # facet_grid(batch~., scales='free') +
  theme_rangeframe() + scale_axis_rangeframe() +
  lims(x=c(-3000,0)) +
  article_text() +
  theme(
    strip.text = element_blank()
  ) +
  labs(x='Log likelihood', y='Count')

ggsave('plots/paper/supp_et/atac_qc_demux_llk_hist1.pdf', width=3.5, height=1.7, units='cm')

ggplot(filter(demux_atac_meta, age==4)) +
  geom_histogram(mapping=aes(llk1), fill=blue, alpha=0.5) +
  geom_histogram(mapping=aes(llk2), fill=redpink, alpha=0.5) +
  scale_y_continuous(breaks=scales::breaks_extended(3)) +
  theme_article() +
  facet_grid(batch~., scales='free') +
  theme_rangeframe() + scale_axis_rangeframe() +
  # lims(x=c(-10000,0)) +
  article_text() +
  theme(
    strip.text = element_blank()
  ) +
  labs(x='Log likelihood', y='Count')
  
ggsave('plots/paper/supp_et/rna_qc_demux_llk_hist1.pdf', width=3.7, height=1.5, units='cm')


rna_atac_meta <- bind_rows('RNA'=demux_rna_meta, 'ATAC'=demux_atac_meta, .id='modality') %>% 
  # filter(!is.na(line)) %>% 
  arrange(age) %>% mutate(
    orig.ident=factor(orig.ident, levels=unique(.$orig.ident)),
    modality=factor(modality, levels=c('RNA', 'ATAC'))
  )

p1 <- ggplot(rna_atac_meta, aes(orig.ident, fill=line)) +
  geom_bar(position='fill') +
  scale_fill_manual(values=line_colors, na.value='grey') +
  facet_grid(~modality, scales='free') +
  scale_y_continuous(breaks=c(0, 1)) +
  theme_rangeframe() + scale_axis_rangeframe() +
  article_text() +
  no_x_text() +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank()
  ) +
  no_legend()
  
  
p2 <- ggplot(rna_atac_meta, aes(orig.ident)) +
  geom_bar(fill='darkgrey', size=0) +
  facet_grid(~modality, scales='free') +
  scale_y_continuous(breaks=c(0, 8000)) +
  theme_rangeframe() + scale_axis_rangeframe() +
  article_text() +
  no_x_text() +
  no_label() +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank()
  )


p2 / p1 + plot_layout(heights=c(1,2)) & no_margin()
ggsave('plots/paper/supp_et/rna_atac_qc_demux_prop.pdf', width=4, height=2, units='cm')  





#### Plots on integration + matching ####

rna_atac_matches <- read_tsv('data/RNA_ATAC/integration/all_matches.tsv')

plot_df <- rna_atac_matches %>% 
  pivot_longer(cols=c(n_cells_atac, n_cells_rna)) %>% 
  mutate(name=factor(name, levels=rev(unique(.$name))))

ggplot(plot_df, aes(value)) +
  geom_histogram(binwidth = 1, color='black', fill='darkgrey', size=0.2) + 
  # scale_y_continuous(breaks=c(0, 20000, 40000)) +
  facet_grid(line~name, scales='free_y') +
  theme_rangeframe() + scale_axis_rangeframe() +
  article_text() +
  theme(
    strip.text = element_blank()
  ) +
  labs(x='Cells per pseudocell', y='Count')

ggsave('plots/paper/supp_et/rna_atac_integration_cell_per_pseudocell_hist.pdf', width=4.7, height=5, units='cm') 



plot_df <- tibble(
  n_matched_rna = length(unique(rna_atac_matches$RNA)),
  n_unmatched_rna = sum(!rna_meta$cell%in%rna_atac_matches$RNA),
  n_matched_atac = length(unique(rna_atac_matches$ATAC)),
  n_unmatched_atac = sum(!atac_meta$cell%in%rna_atac_matches$ATAC)
) %>% pivot_longer(everything()) %>% 
  mutate(
    modality=factor(ifelse(str_detect(name, 'atac'), 'ATAC', 'RNA'), levels=c('RNA', 'ATAC')),
    name=factor(name, levels=rev(unique(.$name)))
  )

ggplot(plot_df, aes(value, name)) +
  geom_bar(stat='identity', color='black', fill='darkgrey', size=0.2) +
  facet_wrap(~modality, scales='free') +
  scale_axis_rangeframe() + theme_rangeframe() +
  scale_x_continuous(breaks=c(0, 25000, 50000), limits=c(0, 60000)) +
  article_text() +
  no_y_text() +
  theme(
    strip.text = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(x='Count')
  
ggsave('plots/paper/supp_et/rna_atac_integration_nmatched_bar.pdf', width=4.4, height=1.5, units='cm') 



#### Check integration: motif/TF expr correlation ####

motif2tf_match <- motif2tf %>% pull(tf, motif)

motif_mat <- GetAssayData(rnatac, assay='chromvar')
gene_expr <- GetAssayData(rnatac, assay='RNA')

tfs_expr <- intersect(motif2tf_match, rownames(gene_expr))
motif2tf_match <- motif2tf_match[match(tfs_expr, motif2tf_match)]
motifs_use <- intersect(names(motif2tf_match), rownames(motif_mat))

tf_mot_mat <- aggregate_matrix(motif_mat[motifs_use, ], groups=motif2tf_match[motifs_use])
tf_expr <- GetAssayData(rnatac, assay='RNA')[unique(motif2tf_match[motifs_use]), ]

motif_tf_cor <- sparse_cov(t(tf_mot_mat), t(tf_expr))$cor

perc_expr <- percent_expression(t(tf_expr)) %>% 
  enframe('TF', 'prc_expr')

motif_tf_cor_df <- motif_tf_cor %>% as_tibble(rownames='tf1') %>% 
  pivot_longer(!tf1, names_to='tf2', values_to='corr') %>% 
  filter(tf1==tf2) %>% 
  select('TF'=tf1, corr) %>% 
  inner_join(perc_expr)

ggplot(motif_tf_cor_df, aes(corr, prc_expr, label=TF)) +
  geom_point() 

goi <- c('FOXG1', 'JUND')
mot_use <- filter(motif2tf, tf%in%goi) %>% pull(motif)
rnatac@active.assay <- 'chromvar'
feature_plot(rnatac, features=mot_use)




#### Gene expression dynamics along pseudotime per line ####
rnatac_meta <- rnatac@meta.data %>% as_tibble(rownames='cell')

goi <- c('POU5F1', 'LIN28A', 'GLI3', 'VIM', 'BCL11B', 'GRIA2')
goi_expr <- GetAssayData(rnatac, assay='RNA')[goi, ] %>% t() %>% 
  as_tibble(rownames='cell') %>% 
  pivot_longer(!cell) %>% 
  inner_join(rnatac_meta) %>% 
  mutate(name=factor(name, levels=goi))



ggplot(goi_expr, aes(pseudotime_ranks, value, color=line)) +
  # geom_point(alpha=0.05, size=0.05) +
  geom_smooth(size=0.3) +
  facet_wrap(~name, scales='free_y', ncol=2) +
  # scale_y_continuous(limits=c(0,2)) +
  theme_void() + no_legend() +
  scale_color_manual(values=line_colors) +
  theme(
    strip.text = element_blank()
  )
ggsave('plots/paper/supp_et/rna_atac_expr_pt_line_smooth.pdf', width=5, height=5, units='cm') 



#### Plot terminal states and VoxHunt plots ####
library(voxhunt)
load_aba_data('/links/groups/treutlein/PUBLIC_DATA/tools/voxhunt_data/')
e13_markers <- structure_markers('E13')
top_genes <- e13_markers %>% 
  group_by(group) %>% 
  top_n(20, auc) %>% pull(gene) %>% unique()


celltype_colors <- c(
  '#E74B3A','#8847A2', '#04597A'
)

celltype_names <- c(
  'ctx_ex', 'ge_in', 'mesen_ex'
)
names(celltype_colors) <- celltype_names

dim_plot(rnatac, group.by='neuron_type') +
  scale_color_manual(values=celltype_colors, na.value='grey') +
  theme_void() + no_legend()
ggsave('plots/paper/supp_et/rna_atac_neuron_types_umap.png', width=3.2*5, height=2*5, units='cm') 

rnatac@active.assay <- 'RNA'
rnatac_neurons <- subset(DietSeurat(rnatac, assays='RNA'), neuron_type%in%celltype_names)

neuron_voxmap <- voxel_map(rnatac_neurons, group_name = 'neuron_type', genes_use=top_genes)
plot_map(neuron_voxmap)


plot_map(neuron_voxmap, view='slice', slices=c(9, 13, 27, 35)) & no_legend() 
ggsave('plots/paper/supp_et/rna_atac_neuron_types_voxmap.pdf', width=2.7, height=6, units='cm') 




#### Feature plots ####
goi <- c('APOE', 'DCN', 'VIM', 'DCX', 'TCF7L2', 'ISL1', 'KRT8', 'COL5A1', 'SOX2', 'NEUROD6', 'LHX9', 'GAD1')
feature_plot(rnatac, features=goi, order=T, pt.size=1, ncol=6) &
  scale_color_gradientn(colors=gyylgnbu()) &
  ggtitle('') &
  no_margin()
ggsave('plots/paper/supp_et/rna_feature_umap.png', width=18*2, height=4.5*2)





#### More feature plots with RNA and ATAC features + graph layout ####

plot_markers <- c('LIN28A', 'GLI3', 'ZIC2', 'FOXG1', 'NFIA', 'NEUROD1', 'NKX2-1', 'DLX5')
rnatac@active.assay <- 'RNA'
pf_ex <- feature_plot(rnatac, features=plot_markers, pt.size=0.1, order=T, ncol=length(plot_markers)) &
  scale_color_gradientn(colors=gyylgnbu(1)) & theme(title = element_blank()) & ggtitle('')

rnatac@active.assay <- 'gene_activity'
pf_act <- feature_plot(rnatac, features=paste0('geneactivity_', plot_markers), pt.size=0.1, order=T, ncol=length(plot_markers)) &
  scale_color_gradientn(colors=gypurd(1)) & theme(title = element_blank()) & ggtitle('')


gene_expr_df <- rnatac_res20_summaries$RNA[, plot_markers] %>% 
  {colnames(.) <- str_replace(colnames(.), '-', '_');.} %>% 
  as_tibble(rownames='name')

perc_expr_df <- rnatac_res20_summaries$perc_expr[, plot_markers] %>% 
  {colnames(.) <- paste0(str_replace(colnames(.), '-', '_'), '_perc');.} %>% 
  as_tibble(rownames='name') 

plot_tree <- lineage_graph %N>% 
  left_join(gene_expr_df) %>% 
  left_join(perc_expr_df)

plots <- map(str_replace(plot_markers, '-', '_'), function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.2, width=0.2) +
    geom_node_point(aes_string(fill=gene, size=paste0(gene, '_perc')), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.3, 2.5)) +
    scale_fill_gradientn(colors=gyylgnbu()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend()
})
p_ex <- wrap_plots(plots, nrow=1)



plot_markers_act <- plot_markers %>% intersect(colnames(rnatac_res20_summaries$gene_activity))
# plot_markers_act <- top_act_de$gene %>% intersect(colnames(cluster_act))
gene_act_df <- rnatac_res20_summaries$gene_activity[, plot_markers_act] %>% 
  {colnames(.) <- str_replace(colnames(.), '-', '_');.} %>% 
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>% 
  left_join(gene_act_df)

plots <- map(str_replace(plot_markers_act, '-', '_'), function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.2, width=0.2) +
    geom_node_point(aes_string(fill=gene), size=2, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=gypurd()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend()
})
p_act <- wrap_plots(plots, nrow=1)

pf_ex / p_ex / pf_act / p_act & no_margin()
ggsave('plots/paper/supp_et/rna_atac_ex_act_feature_umap_graph.png', width=20, height=8)



plot_motif_markers <- c('POU5F1', 'NR2F1', 'OTX2', 'DLX5', 'NFIA')
plot_motifs <- motif2tf %>% 
  filter(tf%in%plot_motif_markers, collection=='CORE') %>% 
  pull(motif, tf)
plot_motifs <- plot_motifs[plot_motif_markers]
plot_markers_mot <- plot_motifs %>% intersect(colnames(rnatac_res20_summaries$chromvar))

rnatac@active.assay <- 'chromvar'
pf_mot <- feature_plot(rnatac, features=paste0('chromvar_', plot_markers_mot), pt.size=0.1, order=T, ncol=length(plot_markers_mot)) &
  scale_color_gradientn(colors=gyylorbr(1)) & theme(title = element_blank()) & ggtitle('')

gene_mot_df <- rnatac_res20_summaries$chromvar[, plot_markers_mot] %>% 
  {colnames(.) <- str_replace(colnames(.), '-', '_');.} %>% 
  as_tibble(rownames='name')

plot_tree <- lineage_graph %N>% 
  left_join(gene_mot_df)

plots <- map(str_replace(plot_markers_mot, '-', '_'), function(gene){
  ggraph(plot_tree, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.2, width=0.2) +
    geom_node_point(aes_string(fill=gene), size=2, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=gyylorbr()) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme_void() + no_legend()
})
p_mot <- wrap_plots(plots, nrow=1)


pf_mot / p_mot & no_margin()
ggsave('plots/paper/supp_et/rna_atac_motif_feature_umap_graph.png', width=12.5, height=4)


print_scale(grad(pals::ocean.tempo, 0.5))
ggsave('plots/paper/supp_et/ocean.tempo.pdf')










