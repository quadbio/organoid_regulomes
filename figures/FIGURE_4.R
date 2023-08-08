source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/perturbator/enrichment.R')
source('~/scripts/grn/models.R')

setwd('~/projects/early/')


pathway_colors <- c('#8DCFBC', '#C5D700', '#B1A8CD', '#E16561', '#7FA1CC', '#EAAB4C', '#EAAB4C')
names(pathway_colors) <- c('Bmp', 'Fgf', 'Hedgehog', 'Notch', 'RA', 'Wnt', 'Wnt-Rspo')

type_colors <- c('#5EB670', '#E6CC07', '#73169E', '#579FFF')
names(type_colors) <- c('Target genes', 'Ligands', 'Receptors', 'Transcription factor')


rnatac %>% feature_plot(features=c('ID1', 'ID2'), order=T)
gli3_45 %>% feature_plot(features=c('ID1', 'BCL11A', 'LHX8', 'NKX2-1', 'DMRTA1'), order=T, split.by='GLI3_KO')



#### Read stuff ####
shh <- read_rds('data/SHH/shh_combined_srt_v1.1preproc.rds')
gli3 <- read_rds('~/projects/CROP_seq/data/GLI3_KO/multiome/DAY18/muo_gli3_all/muo_gli3_all.rds')

gli3_rna <- read_rds('~/projects/CROP_seq/data/GLI3_KO/RNA/GLI3_RNA_all_v2.rds')
# gli3_45 <- read_rds('~/projects/CROP_seq/data/GLI3_KO/RNA/45d/GLI3_RNA_45d.rds')

good_genes <- read_tsv('~/resources/gene_sets/protein_coding_genes_goodlist.txt', col_names=FALSE)$X1
feature_sets <- read_rds('data/gene_sets/RNA_feature_sets.rds')

telen_kowt_de <- read_rds('~/projects/CROP_seq/data/GLI3_KO/multiome/DAY18/muo_gli3_all/diff_expression/telen_ko_wt_de_results.rds')
modules <- read_tsv('data/RNA_ATAC/grn/networks/v5_d100k_cons_cre_noex_modules.tsv')
modules_list <- read_rds('data/RNA_ATAC/grn/networks/tf_peak_gene_modules.rds')
shh_de <- read_tsv('data/SHH/shh_telen_de.tsv')
patterning_genes <- read_tsv('~/resources/gene_sets/human_patterning_genes_v3tfs.tsv')
tfs <- read_tsv('~/resources/DB/animal_tfdb/tf_human.tsv')

gli3_telen_de <- telen_kowt_de$RNA %>%
    filter(group=='TRUE') %>%
    mutate(group='telen_npc')


#### GLI3 DE in 30d and 45d progenitors ####
gli3_rna$state <- ifelse(gli3_rna$seurat_clusters%in%c(34,31,49,21,42,6,0,15,26,30,61,5,1,28,35), 'neuron', 'npc')
gli3_rna$cluster_lineage %>% table(gli3_rna$age)

dim_plot(gli3_rna, group.by=c('state', 'age', 'seurat_clusters'), label=T)

gli3_rna$lineage_state <- paste0(gli3_rna$cluster_lineage, '_', gli3_rna$state)
gli3_rna$lineage_state_age <- paste0(gli3_rna$lineage_state, '_', gli3_rna$age)
gli3_telen <- gli3_rna %>% subset(cluster_lineage%in%c('ge') & age==c('45d'))

gli3_telen_split <- gli3_telen %>% SplitObject(split.by='lineage_state_age')
gli3_kowt_de <- map_dfr(gli3_telen_split, function(x){
    if (length(unique(x$GLI3_status))>1){
        de(x, group='GLI3_status')
    } else {
        return(tibble())
    }
}, .id='state')

gli3_kowt_de <- gli3_kowt_de %>%
    filter(group=='KO') %>%
    select(-group) %>%
    rename('group'='state')

gli3_kowt_de <- gli3_kowt_de %>% filter(avg_exp>0)
# gli3_kowt_de %>% write_tsv('~/projects/CROP_seq/data/GLI3_KO/RNA/ge_lineage_45d_ko_wt_de_results.tsv')
gli3_kowt_de <- read_tsv('~/projects/CROP_seq/data/GLI3_KO/RNA/ge_lineage_45d_ko_wt_de_results.tsv')

combined_de_df <- bind_rows(gli3_telen_de, gli3_kowt_de) %>%
    filter(avg_exp>0)

# combined_de_df %>% write_tsv('~/projects/CROP_seq/data/GLI3_KO/RNA/gli3_telen_ge_ko_wt_de_combined.tsv')
# combined_de_df <- read_tsv('~/projects/CROP_seq/data/GLI3_KO/RNA/gli3_telen_ge_ko_wt_de_combined.tsv')

plot_df <- combined_de_df %>%
    filter(avg_exp>0, padj<1e-4, feature%in%good_genes, abs(fc)>0.1)

vulcano_plot(plot_df, top_only = F)

ggplot(plot_df, aes(group, fc)) +
    geom_line(aes(group=feature), color='grey', alpha=0.5) +
    # geom_point() +
    geom_text(aes(label=feature), size=3)


#### Incorporate patterning genes ####
combined_de_df <- read_tsv('~/projects/CROP_seq/data/GLI3_KO/RNA/gli3_telen_ge_ko_wt_de_combined.tsv')
combined_de_pat <- combined_de_df %>%
    inner_join(patterning_genes, by=c('feature'='symbol')) %>%
    mutate(group=factor(group, levels=c('telen_npc', 'ge_npc_20d', 'ge_npc_30d', 'ge_npc_45d', 'ge_neuron_30d', 'ge_neuron_45d'))) %>%
    filter(group%in%c('telen_npc', 'ge_neuron_45d', 'ge_npc_30d', 'ge_npc_45d')) %>%
    group_by(feature) %>% filter(any(padj<1e-4 & abs(fc)>0.1))
    # group_by(feature) %>% filter(any(padj<1e-4))

combined_df_mat <- combined_de_pat %>%
    select(feature, fc, group) %>%
    distinct() %>%
    pivot_wider(names_from=feature, values_from=fc, values_fill=0) %>%
    column_to_rownames('group') %>% as.matrix()

gene_order <- combined_df_mat %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$labels[.$order]}

early_de <- combined_de_df %>%
    filter(group=='telen_npc', padj<1e-4, abs(fc)>0.1) %>%
    pull(feature)

npc_de <- combined_de_df %>%
    filter(group%in%c('ge_npc_45d'), padj<1e-4, abs(fc)>0.1) %>%
    pull(feature)

neuron_de <- combined_de_df %>%
    filter(group%in%c('ge_neuron_45d'), padj<1e-4, abs(fc)>0.1) %>%
    pull(feature)

pat_de_plot <- combined_de_pat %>%
    mutate(feature=factor(feature, levels=gene_order)) %>%
    mutate(
        fc=ifelse(abs(fc)>1, sign(fc)*1, fc),
        early_de=factor(case_when(
            feature%in%neuron_de & feature%in%npc_de & feature%in%early_de ~ 'early_neuron_npc',
            feature%in%early_de & feature%in%npc_de ~ 'early_npc',
            feature%in%neuron_de & feature%in%npc_de ~ 'neuron_npc',
            feature%in%early_de ~ 'early',
            feature%in%npc_de ~ 'npc',
            feature%in%neuron_de ~ 'neuron'
            ), levels=c('early', 'early_npc', 'npc', 'neuron_npc', 'neuron', 'early_neuron_npc'))
    )

ph <- ggplot(pat_de_plot, aes(group, feature, fill=fc)) +
    geom_tile() +
    scale_fill_gradientn(colors=rev(bigrad(pals::brewer.rdbu, 2)), limits=c(-1,1)) +
    article_text() +
    facet_grid(early_de~., scales='free', space='free') +
    no_y_text() + rotate_x_text(40) +
    no_label() + no_margin() + no_legend() +
    theme(
        panel.border = element_blank(),
        # strip.text = element_blank()
    )
ph
pp <- ggplot(pat_de_plot, aes(feature, fill=pathway)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=pathway_colors, na.value = 'white') +
    facet_grid(early_de~., scales='free', space='free') +
    theme_void() + no_legend() + coord_flip() +
    theme(
        axis.title.y = element_text(size=6, angle=90),
        axis.text.y = element_text(size=5, hjust=1),
        strip.text = element_blank()
    ) +
    labs(x='Gene')

pt <- ggplot(pat_de_plot, aes(feature, fill=type)) +
    facet_grid(early_de~., scales='free', space='free') +
    geom_bar(position='fill') +
    scale_fill_manual(values=type_colors, na.value = 'white') +
    theme_void() + no_legend() + coord_flip() +
    theme(
        strip.text = element_blank()
    )


(pp | pt | ph) + plot_layout(widths=c(1,1,10))
ggsave('plots/paper/fig4/fig4_gli3_de_pat_heatmap2.pdf', width=3.5, height=14, unit='cm')


#### SHH UMAP ####
shh %>% dim_plot(reduction='wnnumap', label=T)
shh_use <- shh %>% subset(seurat_clusters%in%c(7,8,0,13,1,0,16,10,11,2,3,9,5,12))

shh_use <- shh_use %>% 
    FindMultiModalNeighbors(reduction.list = c('pca', 'lsi'), dims.list=list(1:20,2:20)) 
shh_use <- shh_use %>% 
    RunUMAP(nn.name = 'weighted.nn')
p1 <- shh_use %>% dim_plot(reduction='umap', label=T)
p2 <- shh_use %>% feature_plot(reduction='umap', features=c('WLS', 'FOXG1'))

p1 | p2

cluster_colors <- c(
    '#883D79', '#E0BAD8', '#652D59', '#6B4283',
    '#A8E4FC', '#05668d', '#2C86BD',
    '#3A3D54', '#d0d3d4', '#686A6A', '#DDDDDD', '#BBBBBB', '#555555', '#777777'
)

names(cluster_colors) <- c(
    7,8,0,9,
    10,11,2,
    16,13,1,3,9,5,12
)


cluster_umap_means <- Reductions(shh_use, 'umap')@cell.embeddings %>%
    aggregate_matrix(as.character(shh_use$seurat_clusters)) %>% as_tibble(rownames='seurat_clusters')
plot_df <- shh_use@meta.data %>%
    as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(shh_use[['umap']]@cell.embeddings, rownames='cell')) %>% 
    filter(seurat_clusters%in%names(cluster_colors))
ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=seurat_clusters)) +
    geom_point(size=1, shape=21, stroke=0.1) +
    scale_fill_manual(values=cluster_colors) +
    theme_void() + no_legend()

ggsave('plots/paper/fig4/fig4_shh_clusters_umap.png', width=6, height=5, unit='cm')


ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=orig.ident)) +
    geom_point(size=1, shape=21, stroke=0.05) +
    scale_fill_manual(values=c('grey', 'black')) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/fig4_shh_condition_umap.png', width=6, height=5, unit='cm')

shh_use %>% feature_plot(reduction='umap', features=c('WLS', 'FOXG1', 'POU5F1', 'GLI3', 'FGF8', 'PAX3'), order=T) &
    scale_color_gradientn(colors=gyylgnbu())
ggsave('plots/paper/fig4/fig4_shh_features_umap.png', width=12, height=18, unit='cm')



#### Plot comparison between GLI3 and SHH ####
gli3_module <- modules %>% filter(tf=='GLI3')

shh_telen_compare <- shh_de %>%
    full_join(gli3_telen_de, by=c('feature')) %>%
    left_join(gli3_module, by=c('feature'='target')) %>%
    mutate(grn_dir=sign(estimate)) %>%
    mutate(sig_status=case_when(
        padj.x < 1e-4 & padj.y < 1e-4 ~ 'both',
        padj.x < 1e-4 & padj.y >= 1e-4 ~ 'SHH',
        padj.x >= 1e-4 & padj.y < 1e-4 ~ 'GLI3',
        T ~ 'none'
    ))


cor(shh_telen_compare$fc.x, shh_telen_compare$fc.y)

plot_df <- shh_telen_compare %>%
    filter(padj.x < 0.1 | padj.y < 0.1)

sig_colors <- c('GLI3'='#668795', 'SHH'='#B49AC4', 'both'='black')

ggplot(plot_df, aes(fc.x, fc.y, fill=factor(sig_status))) +
    geom_hline(yintercept = 0, color='darkgrey', size=0.2) +
    geom_vline(xintercept = 0, color='darkgrey', size=0.2) +
    geom_point(data=filter(plot_df, sig_status=='none'), fill='grey', alpha=0.3, shape=21, stroke=0.05, size=0.5) +
    geom_point(data=filter(plot_df, sig_status!='none'), shape=21, stroke=0.05, size=0.5) +
    scale_fill_manual(values=sig_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    # geom_text_repel(data=filter(shh_telen_compare, abs(fc.y)>0.3 | abs(fc.x)>0.3), mapping=aes(label=feature), max.overlaps=99999, size=2) +
    no_legend() +
    labs(x='Fold change (SHH vs control)', y='Fold change (GLI3 KO vs WT)')
ggsave('plots/paper/fig4/SHH_vs_GLI3_de.pdf', width=4.5, height=4, unit='cm')



ggplot(plot_df, aes(fc.x, fc.y, fill=factor(sig_status))) +
    geom_hline(yintercept = 0, color='darkgrey', size=0.2) +
    geom_vline(xintercept = 0, color='darkgrey', size=0.2) +
    geom_point(data=filter(plot_df, sig_status=='none'), fill='grey', alpha=0.3, shape=21, stroke=0.05, size=0.5) +
    geom_point(data=filter(plot_df, sig_status!='none'), shape=21, stroke=0.05, size=0.5) +
    scale_fill_manual(values=sig_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    geom_text_repel(data=filter(shh_telen_compare, abs(fc.y)>0.3 | abs(fc.x)>0.3), mapping=aes(label=feature), max.overlaps=99999, size=2) +
    no_legend() +
    labs(x='Fold change (SHH vs control)', y='Fold change (GLI3 KO vs WT)')




#### GO enrichment of sig sets ####
shh_gli3_sig <- shh_telen_compare %>%
    mutate(
        gli3_sig = (padj.y<1e-4) & (abs(fc.y)>0.1),
        shh_sig = (padj.x<1e-4) & (abs(fc.x)>0.1)
    )

perc_expr <- rowMeans(shh[['RNA']]@data > 0)
expr_genes <- names(perc_expr[perc_expr>0.05])

all_features <- bitr(expr_genes, fromType = 'SYMBOL', toType = c('ENSEMBL', 'ENTREZID'), OrgDb = org.Hs.eg.db) %>%
    as_tibble()

shh_deg <- filter(shh_gli3_sig, shh_sig & !gli3_sig, feature%in%good_genes)$feature
shh_deg_tfs <- shh_deg %>% intersect(modules$tf)

gli3_deg <- filter(shh_gli3_sig, !shh_sig & gli3_sig, feature%in%good_genes)$feature
common_deg <- filter(shh_gli3_sig, shh_sig & gli3_sig, sign(fc.x)==sign(fc.y), feature%in%good_genes)$feature
shh_opposite_deg <- filter(shh_gli3_sig, shh_sig & gli3_sig, sign(fc.x)!=sign(fc.y), feature%in%good_genes)$feature


shh_gene_ids <- filter(all_features, SYMBOL%in%shh_deg)
shh_ego <- enrichGO(
    gene = shh_gene_ids$ENTREZID,
    universe = unique(all_features$ENTREZID),
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    # qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

gli3_gene_ids <- filter(all_features, SYMBOL%in%gli3_deg)
gli3_ego <- enrichGO(
    gene = gli3_gene_ids$ENTREZID,
    universe = unique(all_features$ENTREZID),
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    # qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

common_gene_ids <- filter(all_features, SYMBOL%in%common_deg)
common_ego <- enrichGO(
    gene = common_gene_ids$ENTREZID,
    universe = unique(all_features$ENTREZID),
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    # qvalueCutoff = 0.05,
    pAdjustMethod = 'fdr',
    ont = 'ALL',
    readable = T
)

all_ego <- bind_rows(
    'SHH' = as_tibble(shh_ego@result),
    'GLI3' = as_tibble(gli3_ego@result),
    'both' = as_tibble(common_ego@result),
    .id='status'
) %>%
    mutate(
        n_in=as.numeric(str_replace(GeneRatio, '(\\d+)/\\d+', '\\1')),
        n_out=as.numeric(str_replace(GeneRatio, '\\d+/(\\d+)', '\\1')),
        bg_in=as.numeric(str_replace(BgRatio, '(\\d+)/\\d+', '\\1')),
        bg_out=as.numeric(str_replace(BgRatio, '\\d+/(\\d+)', '\\1')),
        odds_ratio = (n_in/n_out) / (bg_in/bg_out),
        log_or = log2(odds_ratio)
    )

all_ego  %>% write_tsv('data/SHH/shh_gli3_de_GO_enrich.tsv')

plot_df <- all_ego %>%
    filter(qvalue<1e-2, ONTOLOGY=='BP') %>%
    group_by(status) %>%
    top_n(30, log_or)


ggplot(plot_df, aes(log_or, Description)) +
    geom_bar(stat='identity') +
    facet_grid(status~., scales='free')





#### Check lineage prob correlations of GLI3 downstream modules / all DE TFs ####
rnatac <- read_rds('/local2/USERS/jfleck/projects/early/data/RNA_ATAC_pseudocells_v2.2links_srt.rds')
lineage_graph <- read_rds('~/projects/early/data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')
hc_meta <- read_tsv('~/projects/early/data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv') %>%
    mutate(cluster=as.character(cluster))
cluster_meta <- as_tibble(lineage_graph)
gene_meta <- read_tsv('~/projects/early/data/gene_sets/gene_scores.tsv')
spath_graph_agree <- read_tsv('~/projects/CROP_seq/data/GLI3_KO/multiome/DAY18/muo_gli3_all/grn/gli3_consistent_de_spaths.tsv')
branch_modules <- read_tsv('data/RNA_ATAC/grn/networks/v5_branchpoint_modules.tsv')

hc_meta_use <- cluster_meta %>%
    filter(lineage!='early')

ctx_rna_corr <- cor(hc_meta_use$to_ctx_ranks, as.matrix(rnatac@misc$cluster_summaries$RNA[hc_meta_use$name, ]), method = 'pearson')
ge_rna_corr <- cor(hc_meta_use$to_ge_ranks, as.matrix(rnatac@misc$cluster_summaries$RNA[hc_meta_use$name, ]), method = 'pearson')
nt_rna_corr <- cor(hc_meta_use$to_nt_ranks, as.matrix(rnatac@misc$cluster_summaries$RNA[hc_meta_use$name, ]), method = 'pearson')

ctx_rna_corr_scale <- scale(ctx_rna_corr[1,])
ge_rna_corr_scale <- scale(ge_rna_corr[1,])
nt_rna_corr_scale <- scale(nt_rna_corr[1,])

ctx_rna_score <- ctx_rna_corr %>% {.[.<0]<-0;.}
ge_rna_score <- ge_rna_corr %>% {.[.<0]<-0;.}
nt_rna_score <- nt_rna_corr %>% {.[.<0]<-0;.}


ctx_rna_corr_scale %>% hist()
ge_rna_corr_scale %>% hist()
nt_rna_corr_scale %>% hist()

lineage_probs <- do.call(rbind, list(ctx_rna_score, ge_rna_score, nt_rna_score)) %>% t() %>%
    as_tibble(rownames='gene') %>% {colnames(.)[-1] <- c('ctx_score', 'ge_score', 'nt_score'); .} %>%
    filter(!is.na(ge_score))

lineage_probs <- do.call(rbind, list(ctx_rna_corr, ge_rna_corr, nt_rna_corr)) %>% t() %>%
    as_tibble(rownames='gene') %>% {colnames(.)[-1] <- c('ctx_score', 'ge_score', 'nt_score'); .} %>%
    filter(!is.na(ge_score))

lineage_scores <- do.call(cbind, list(ctx_rna_corr_scale, ge_rna_corr_scale, nt_rna_corr_scale)) %>%
    as_tibble(rownames='gene') %>% {colnames(.)[-1] <- c('ctx_score', 'ge_score', 'nt_score'); .} %>%
    filter(!is.na(ge_score))

gli3_down_tfs <- filter(spath_graph_agree, from_node=='GLI3')$to_node %>% intersect(modules$tf)
gli3_down_tf_modules <- filter(modules, tf%in%gli3_down_tfs)

telen_de_df <- combined_de_df %>%
    filter(group=='telen_npc', padj<1e-4, abs(fc)>0.1, feature%in%modules$tf)

gli3_de_tf_modules <- filter(branch_modules) %>%
    inner_join(telen_de_df, by=c('tf'='feature'))

lineage_path_cor <- lineage_scores %>%
    # left_join(patterning_genes, by=c('gene'='symbol')) %>%
    inner_join(gli3_de_tf_modules, by=c('gene'='target')) %>%
    mutate(dir=sign(estimate))


plot_df <- lineage_path_cor %>%
    filter(dir==1) %>%
    group_by(tf) %>%
    mutate(corr_mean=median(ctx_score)) %>%
    arrange(desc(corr_mean)) %>%
    mutate(tf=factor(tf, levels=unique(.$tf)))

p1 <- ggplot(plot_df, aes(tf, ctx_score, fill=factor(sign(fc)))) +
    geom_hline(yintercept = 0, color='grey', size=0.3) +
    # geom_quasirandom(color='darkgrey', alpha=0.8, size=0.1) +
    geom_violin(fill='darkgrey', alpha=0.8, size=0.1) +
    geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, size=0.2, fill=lineage_colors['ctx']) +
    article_text() +
    theme_rangeframe() + scale_axis_rangeframe() +
    rotate_x_text(40) +
    # coord_flip() +
    labs(x='TF module', y='Correlation with dorsal telen.')


plot_df <- lineage_path_cor %>%
    filter(dir==1) %>%
    group_by(tf) %>%
    mutate(corr_mean=median(ge_score)) %>%
    arrange(desc(corr_mean)) %>%
    mutate(tf=factor(tf, levels=unique(.$tf)))

p2 <- ggplot(plot_df, aes(tf, ge_score, fill=factor(sign(fc)))) +
    geom_hline(yintercept = 0, color='grey', size=0.3) +
    # geom_quasirandom(color='darkgrey', alpha=0.8, size=0.1) +
    geom_violin(fill='darkgrey', alpha=0.8, size=0.1) +
    geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, size=0.2, fill=lineage_colors['ge']) +
    article_text() +
    theme_rangeframe() + scale_axis_rangeframe() +
    rotate_x_text(40) +
    # coord_flip() +
    labs(x='TF module', y='Correlation with dorsal telen.')

p1 / p2

ggsave('plots/paper/fig4/fig4_gli3_de_tf_ctx_corr_boxplot.pdf', width=6.8, height=5.4, unit='cm')



#### Summarize to TF -> activation/repression score ####
lineage_score_summary <- lineage_path_cor %>%
    group_by(tf, dir) %>%
    summarize(
        ge_mean=mean(ge_score),
        ctx_mean=mean(ctx_score),
        ge_sd=sd(ge_score),
        ctx_sd=sd(ctx_score),
        fc=fc[1],
        fc=ifelse(abs(fc)>1, sign(fc)*1, fc),
        ge_p=pnorm(ge_mean),
        ctx_p=pnorm(ctx_mean)
    ) %>%
    filter(dir==1)

plot_df <- lineage_score_summary %>%
    arrange(desc(ctx_mean)) %>%
    mutate(tf=factor(tf, levels=unique(.$tf)))

p1 <- ggplot(plot_df, aes(tf, ctx_mean, fill=ctx_mean)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    geom_point(shape=21, size=0.6, stroke=0.2) +
    scale_fill_gradientn(colors=bigrad(pals::brewer.reds, 2)) +
    scale_y_continuous(breaks=seq(-1.5,1,0.5)) +
    labs(y='Activation score')


plot_df <- lineage_score_summary %>%
    arrange(desc(ge_mean)) %>%
    mutate(tf=factor(tf, levels=unique(.$tf)))

p2 <- ggplot(plot_df, aes(tf, ge_mean, fill=ge_mean)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    geom_point(shape=21, size=0.6, stroke=0.2) +
    scale_fill_gradientn(colors=bigrad(pals::brewer.purples, 2)) +
    theme(
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
    )


(p1 | p2) & no_legend() & no_margin()
ggsave('plots/paper/fig4/fig4_tf_act_score_lolli.pdf', width=5.5, height=2, units='cm')


plot_df <- lineage_score_summary %>%
    arrange(desc(ctx_mean)) %>%
    mutate(tf=factor(tf, levels=unique(.$tf)))

p1 <- ggplot(plot_df, aes(tf, ctx_mean, fill=ctx_mean)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    rotate_x_text(90) +
    geom_point(shape=21, size=0.6, stroke=0.2) +
    scale_fill_gradientn(colors=bigrad(pals::brewer.reds, 2)) +
    scale_y_continuous(breaks=seq(-1.5,1,0.5)) +
    labs(y='Activation score')


plot_df <- lineage_score_summary %>%
    arrange(desc(ge_mean)) %>%
    mutate(tf=factor(tf, levels=unique(.$tf)))


p2 <- ggplot(plot_df, aes(tf, ge_mean, fill=ge_mean)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    rotate_x_text(90) +
    geom_point(shape=21, size=0.6, stroke=0.2) +
    scale_fill_gradientn(colors=bigrad(pals::brewer.purples, 2)) +
    theme(
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
    )

p1 + p2


#### GLI3 CRE DA ####
shh_gli3 <- read_rds('data/SHH/shh_gli3_muo_srt_v1.1links.rds')
shh_gli3_da <- read_tsv('data/SHH/shh_gli3_telen_da.tsv')

links_df <- Links(shh_gli3) %>%
    as_tibble()

da_links <- shh_gli3_da %>%
    mutate(padj=p.adjust(pval, method='fdr')) %>%
    # filter(padj<1e-2) %>%
    inner_join(links_df, by=c('feature'='peak'))

gli3_da <- filter(da_links, dataset=='GLI3')


#### Check overlaps with H3K27 marks ####
H3K27 <- read_rds('data/CT_H3K27ac_d19_peaks.rds')

ov_peaks <- gli3_da$feature[subjectHits(findOverlaps(H3K27, StringToGRanges(gli3_da$feature)))]


#### Check overlaps with human specific peaks ####
human_chimp_peaks <- read.table('/links/groups/treutlein/USERS/zhisong_he/Work/brain_organoid_variability/analysis_chimp/organoids/analysis/ATACseq/data.peak_DA-DE_info_merged-peaks_transition-as-neuron_with-read-subsample.tsv') %>%
    as_tibble(rownames='peak_name') %>% mutate(peak=paste(chr, start, end, sep='-'))

hg19_to_hg38_chain <- import.chain('~/resources/hg19ToHg38.over.chain')
hg19_peaks <- StringToGRanges(human_chimp_peaks$peak)
hg19_peaks$peak_name <- human_chimp_peaks$peak
hg38_peaks <- liftOver(hg19_peaks, hg19_to_hg38_chain) %>% unlist()

human_chimp_peaks$peak_hg38 <- GRangesToString(hg38_peaks)[match(human_chimp_peaks$peak, hg38_peaks$peak_name)]


# Get human-high peaks
human_npc_da <- human_chimp_peaks %>%
    filter(DA_sign==1, !is.na(peak_hg38)) %>%
    pull(peak_hg38) %>% StringToGRanges()

spec_ov_peaks <- gli3_da$feature[subjectHits(findOverlaps(human_npc_da, StringToGRanges(gli3_da$feature)))]


#### Join with marks ####
gli3_da <- gli3_da %>%
    mutate(has_H3K27ac=feature%in%ov_peaks, is_human_specific=feature%in%spec_ov_peaks)

spec_da <- gli3_da %>%
    filter(is_human_specific) %>%
    group_by(gene) %>%
    filter(padj==min(padj))

mark_da <- gli3_da %>%
    filter(!gene%in%spec_da$gene) %>%
    filter(has_H3K27ac) %>%
    group_by(gene) %>%
    filter(padj==min(padj))

nomark_da <- gli3_da %>%
    filter(!gene%in%mark_da$gene & !gene%in%spec_da$gene) %>%
    group_by(gene) %>%
    filter(padj==min(padj))

combined_da_pat <- bind_rows(mark_da, nomark_da, spec_da) %>%
    filter(padj<1e-4) %>%
    right_join(pat_de_plot, by=c('gene'='feature'))

pat_cre_plot <- combined_da_pat %>%
    filter(gene%in%gene_order) %>%
    mutate(gene=factor(gene, levels=gene_order)) %>%
    distinct(gene, log_dr, early_de, has_H3K27ac, is_human_specific) %>%
    mutate(log_dr=ifelse(abs(log_dr)>3, sign(log_dr)*3, log_dr))


ph <- ggplot(pat_cre_plot, aes('', gene, fill=log_dr)) +
    geom_tile() +
    scale_fill_gradientn(colors=rev(bigrad(pals::brewer.puor, 2)), limits=c(-4,4), na.value = 'white') +
    theme_void() + no_legend() +
    facet_grid(early_de~., scales='free', space='free') +
    theme(
        panel.border = element_blank(),
        strip.text = element_blank()
    )

pm <- ggplot(pat_cre_plot, aes('', gene, fill=has_H3K27ac)) +
    geom_tile() +
    scale_fill_manual(values=c('TRUE'='#424949', 'FALSE'='grey'), na.value = 'white') +
    theme_void() + no_legend() +
    facet_grid(early_de~., scales='free', space='free') +
    theme(
        panel.border = element_blank(),
        strip.text = element_blank()
    )

pe <- ggplot(pat_cre_plot, aes('', gene, fill=is_human_specific)) +
    geom_tile() +
    scale_fill_manual(values=c('TRUE'='#424949', 'FALSE'='grey'), na.value = 'white') +
    theme_void() + no_legend() +
    facet_grid(early_de~., scales='free', space='free') +
    theme(
        panel.border = element_blank(),
        strip.text = element_blank()
    )

ph | pm | pe
ggsave('plots/paper/fig4/fig4_gli3_da_pat_heatmap.pdf', width=1, height=10, unit='cm')



print_scale(rev(bigrad(pals::brewer.puor, 2)))
ggsave('plots/paper/fig4/fig4_puor.pdf', width=5, height=5, unit='cm')


print_scale(rev(bigrad(pals::brewer.rdbu, 2)))
ggsave('plots/paper/fig4/fig4_rdbu.pdf', width=5, height=5, unit='cm')





#### GLI3-centered GRN with all RNA DA and patterning labels ####
combined_de_df <- read_tsv('~/projects/CROP_seq/data/GLI3_KO/RNA/gli3_telen_ge_ko_wt_de_combined.tsv')
module_nodes <- union(modules$target, modules$tf)
module_graph <- modules %>% as_tbl_graph()

ge_npc_deg <- combined_de_df %>%
    filter(group%in%c('ge_npc_45d', 'ge_neuron_45d')) %>%
    filter(avg_exp>0, padj<1e-4, feature%in%good_genes, abs(fc)>0.1, feature%in%module_nodes)

telen_deg <- combined_de_df %>%
    filter(group%in%c('telen_npc')) %>%
    filter(avg_exp>0, padj<1e-4, feature%in%good_genes, abs(fc)>0.1, feature%in%module_nodes)

combined_de_genes <- combined_de_df %>%
    filter(avg_exp>0, padj<1e-4, feature%in%good_genes, abs(fc)>0.1, feature%in%module_nodes) %>%
    pull(feature) %>% unique() %>% sort()


#### Get new spaths to ventral NPC DEG ####
deg_use <- ge_npc_deg
deg_use <- telen_deg

de_genes <- deg_use %>%
    pull(feature) %>% unique() %>% sort()

spaths <- all_shortest_paths(module_graph, 'GLI3', de_genes, mode='out')$res
spath_list <- map_par(spaths, function(p){
    edg <- names(p)
    edg_graph <- module_graph %>%
        filter(name%in%edg) %>%
        convert(to_shortest_path, from=which(.N()$name==edg[1]), to=which(.N()$name==edg[length(edg)])) %E>%
        mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
        as_tibble()
    edg_dir <- edg_graph %>% pull(estimate) %>% sign() %>% prod()
    edg_p <- edg_graph %>% pull(log_padj) %>% mean()
    return(
        list(
            path = tibble(
                start_node = edg[1],
                end_node = edg[length(edg)],
                dir = edg_dir,
                path = paste(edg, collapse=';'),
                order = length(edg)-1,
                mean_padj = edg_p
            ),
            graph = mutate(edg_graph, path=paste(edg, collapse=';'), end_node=edg[length(edg)], comb_dir=edg_dir)
        )
    )
})

spath_dir <- map_dfr(spath_list, function(x) x$path) %>%
    mutate(path_genes=str_split(path, ';'))
spath_graph <- map_dfr(spath_list, function(x) x$graph)


#### Filter DEG only ####
sig_paths <- spath_dir %>% pull(path_genes, path) %>% {.[map_lgl(., function(x) all(x%in%combined_de_genes))]} %>% names()

de_dir <- deg_use %>%
    mutate(dir=-sign(fc)) %>%
    pull(dir, feature)

order <- 3
spath_use <- filter(spath_dir, end_node%in%de_genes, order<=order, path%in%sig_paths)
spath_agree <- spath_use[de_dir[spath_use$end_node] == spath_use$dir, ]

o1_agree <- spath_agree %>%
    filter(order==1)

o2_agree <- spath_agree %>%
    filter(str_detect(path, paste(o1_agree$path, collapse='|')), order==2)

o3_agree <- spath_agree %>%
    filter(str_detect(path, paste(o1_agree$path, collapse='|')), order==3)

module_agree_pruned <- bind_rows(o1_agree, o2_agree) %>%
    select(start_node, end_node, everything()) %>%
    group_by(end_node) %>% filter(mean_padj==max(mean_padj))
    # filter(end_node%in%tfs$symbol | end_node%in%patterning_genes$symbol)

spath_graph_agree <- spath_graph %>%
    filter(path%in%module_agree_pruned$path) %>%
    select(from_node, to_node, end_node, comb_dir) %>% distinct()

patterning_genes_dist <- patterning_genes %>%
    group_by(symbol) %>% filter(row_number()==1)

module_agree_graph <- module_graph %E>%
    mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
    as_tibble() %>% distinct() %>%
    inner_join(spath_graph_agree) %>%
    select(from_node, to_node, everything(), -from, -to) %>% arrange(comb_dir) %>% as_tbl_graph() %N>%
    left_join(patterning_genes_dist, by=c('name'='symbol'))


ggraph(module_agree_graph, layout='tree') +
    geom_edge_diagonal(aes(color=sign(estimate)), alpha=0.5, width=0.2) +
    geom_node_point(aes(fill=pathway), size=1, shape=21, stroke=0.2) +
    geom_node_label(aes(label=name, fill=pathway), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05) +
    scale_edge_color_gradientn(colors=c('#B94D60', '#70AD8C')) +
    scale_fill_manual(values=pathway_colors, na.value = 'white') +
    scale_x_continuous(expand=c(0.1, 0)) +
    scale_y_continuous(expand=c(0.1, 0), trans = 'reverse') +
    article_text() + theme_void() + no_legend() +
    coord_flip()
ggsave('plots/paper/fig4/fig4_gli3_grn_agree_telen_tree.png', width=5, height=25, units='cm')


ggraph(module_agree_graph, layout='tree', circular=T) +
    geom_edge_diagonal(aes(color=sign(estimate)), alpha=0.5, width=0.2) +
    geom_node_point(aes(fill=pathway), size=1, shape=21, stroke=0.2, fill='grey') +
    geom_node_label(aes(label=name, filter=name%in%tfs$symbol), size=4/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05) +
    scale_edge_color_gradientn(colors=c('#B94D60', '#70AD8C')) +
    scale_fill_manual(values=pathway_colors, na.value = 'white') +
    scale_x_continuous(expand=c(0.1, 0)) +
    scale_y_continuous(expand=c(0.1, 0), trans = 'reverse') +
    article_text() + theme_void() + no_legend() +
    coord_flip()
# ggsave('plots/paper/fig4/fig4_gli3_grn_agree_ge_npc_circular.png', width=2.5, height=2.5, units='cm')
# ggsave('plots/paper/fig4/fig4_gli3_grn_agree_ge_npc_circular.pdf', width=2.5, height=2.5, units='cm')
#
# ggsave('plots/paper/fig4/fig4_gli3_grn_agree_telen_circular.png', width=5, height=5, units='cm')
# ggsave('plots/paper/fig4/fig4_gli3_grn_agree_telen_circular.pdf', width=5, height=5, units='cm')
#
ggsave('plots/paper/fig4/fig4_gli3_grn_agree_ge_circular.png', width=5, height=5, units='cm')
ggsave('plots/paper/fig4/fig4_gli3_grn_agree_ge_circular.pdf', width=5, height=5, units='cm')



#### GLI3 binding at 3 weeks vs DE ####
library(rtracklayer)
gli3_snap <- rtracklayer::import('/links/groups/treutlein/DATA/sequencing/20210604_P1587_FIDES_CUT_TAG_bulk_CT_div_Gli3/bamcoverage/210511_29_FZ_bCT_Gli3_native_4C_23d_bin200.bw')
# gli3_snap <- rtracklayer::import('/local1/DATA/sequencing/20220211_P1846_CUT_TAG_bCT/bamcoverage/220202_16_FZ_bCT_Gli3_org1_A4_19d_noSHH_bin200.bw')
gli3_peaks <- rtracklayer::import('/links/groups/treutlein/DATA/sequencing/-/narrow_peaks/narrow_peaks/210511_29_Gli3_native_narrow_peaks_peaks.narrowPeak')

gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')
gene_annot <- CollapseToLongestTranscript(gene_annot)

gene_annot_ext <- Extend(gene_annot, upstream = 2000)

gene_gli3_overlaps <- findOverlaps(gene_annot_ext, gli3_snap)
gene_gli3_score <- gli3_snap$score[subjectHits(gene_gli3_overlaps)] %>%
    as.matrix() %>%
    Pando::aggregate_matrix(gene_annot_ext$gene_name[queryHits(gene_gli3_overlaps)], fun = 'sum')

gene_score <- gene_annot %>%
    as_tibble() %>%
    inner_join(as_tibble(gene_gli3_score, rownames='gene_name')) %>%
    rename('ct_score'='V1') %>%
    mutate(log_score=log10(ct_score))

telen_cre_da <- gli3_da %>%
    filter(padj<1e-4) %>%
    group_by(gene) %>%
    filter(padj==min(padj))

telen_de <- combined_de_df %>%
    filter(group=='telen_npc')

plot_df <- telen_de %>%
    left_join(telen_cre_da, by=c('feature'='gene')) %>%
    left_join(gene_score, by=c('feature'='gene_name')) %>%
    filter(feature%in%good_genes) %>%
    mutate(
        score_scale=as.numeric(scale(ct_score)),
        score_p=pnorm(score_scale, lower.tail = F),
        score_clip=ifelse(ct_score>4000, 4000, ct_score),
        score_scale_clip=ifelse(score_scale>10, 10, score_scale),
        score_norm=(ct_score/width.y),
        score_rank=rank(ct_score)/max(rank(ct_score)),
        high=padj.x<1e-4&abs(fc)>0.1&score_rank>0.9
    )

ggplot(plot_df, aes(fc, score_clip, label=feature, fill=!is.na(log_dr), color=!is.na(log_dr))) +
    geom_point(data=filter(plot_df, !high), shape=16, size=0.2, alpha=0.2) +
    geom_point(data=filter(plot_df, high), shape=21, stroke=0.1, size=0.6, alpha=1, color='black') +
    geom_text_repel(
        data=filter(
            plot_df,
            feature%in%c('HES4', 'HES5', 'HES1', 'PAX6', 'CREB5', 'DMRTA1', 'GSE1', 'PTCH1', 'AUTS2')),
        alpha=1,
        color='black',
        size=5/ggplot2::.pt,
        max.overlaps=99999
        ) +
    scale_fill_manual(values=c('grey', 'black')) +
    scale_color_manual(values=c('grey', 'black')) +
    scale_alpha_manual(values=c(0.1, 1)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_legend() +
    labs(y='GLI3 binding score', x='Log fold change')
ggsave('plots/paper/fig4/fig4_gli3_cnt_vs_fc_scatter.pdf', width=4, height=3.2, units='cm')




ggplot(plot_df, aes(fc, score_clip, label=feature, fill=!is.na(log_dr), color=!is.na(log_dr))) +
    geom_point(data=filter(plot_df, !high), shape=16, size=0.2, alpha=0.2) +
    geom_point(data=filter(plot_df, high), shape=21, stroke=0.1, size=0.6, alpha=1, color='black') +
    geom_text(
        data=filter(plot_df, high),
        alpha=1,
        color='black',
        size=6/ggplot2::.pt,
    ) +
    scale_fill_manual(values=c('grey', 'black')) +
    scale_color_manual(values=c('grey', 'black')) +
    scale_alpha_manual(values=c(0.1, 1)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_legend() +
    labs(y='GLI3 binding score', x='Log fold change')



