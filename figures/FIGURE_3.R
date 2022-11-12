source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/crop.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')
source('~/scripts/perturbator/ko_inference.R')
source('~/scripts/perturbator/plots.R')

# setwd('~/projects/early/')
setwd('~/projects/CROP_seq/')

library(ggraph)
library(tidygraph)

registerDoParallel(cores=50)

#### Load data ####
crop <- read_rds('data/CROP_3/crop_3_all/crop_3_withguides_v2.rds')
# crop %>% write_rds('data/CROP_3/crop_3_all/crop_3_withguides_v2.rds')

rnatac <- read_rds('/local2/USERS/jfleck/projects/early/data/RNA_ATAC_pseudocells_v2.2links_srt.rds')
lineage_graph <- read_rds('~/projects/early/data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')
lineage_tree <- read_rds('~/projects/early/data/RNA_ATAC/lineages/cellrank_lineage_tree.rds')
hc_meta <- read_tsv('~/projects/early/data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv') %>%
    mutate(cluster=as.character(cluster))
cluster_meta <- as_tibble(lineage_graph)

gene_meta <- read_tsv('~/projects/early/data/gene_sets/gene_scores.tsv')

crop_all <- read_rds('data/CROP_3/crop_3_all/crop_3_all.rds')

crop_guide_counts <- read_tsv('data/CROP_3/crop_3_all/guides_new/crop3_guide_counts.tsv') %>%
    as_matrix()

crop <-  add_guide_assay(crop, crop_guide_counts)
crop_all <-  add_guide_assay(crop_all, crop_guide_counts)

crop[['umap']] <- CreateDimReducObject(as.matrix(crop@meta.data[, c('UMAP1', 'UMAP2')]), key='UMAP_')


crop$lineage <- case_when(
    str_detect(crop$celltype, 'ctx') ~ 'ctx',
    str_detect(crop$celltype, 'ge') ~ 'ge',
    str_detect(crop$celltype, 'mesen') ~ 'nt'
)

tf_graph <- read_rds('~/projects/early/data/RNA_ATAC/grn/tf_graph.rds')
tf_graph_meta <- tf_graph %N>% as_tibble()

modules <- read_tsv('~/projects/early/data/RNA_ATAC/grn/networks/v5_d100k_cons_cre_noex_modules.tsv')
module_meta <- modules %>% 
    group_by(tf) %>% 
    summarize(
        ngenes=ngenes[1],
        npeaks=length(unlist(str_split(peaks, ';')))
    ) %>% rename('gene'='tf')

## Peak matches
peak_subsets <- read_rds('~/projects/early/data/RNA_ATAC/RNA_ATAC_peak_subsets_v2.rds')
grn_regions <- peak_subsets$cons_cre_noex
grnreg2peaks <- rownames(rnatac@assays$peaks@data)[grn_regions$peaks]
names(grnreg2peaks) <- GRangesToString(grn_regions$ranges)


#### Annotation and feature plots ####
crop <- CellCycleScoring(crop, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)

meta <- crop@meta.data %>%
  as_tibble(rownames='cell') %>%
  arrange(factor(org_group))
ggplot(meta, aes(-UMAP1, -UMAP2, fill=org_group)) +
  geom_point(size=1, shape=21, stroke=0.1) +
  theme_void() + no_legend() +
  scale_fill_manual(values=c('#eceff1', '#b0bec5', '#546e7a', '#263238'))
ggsave('plots/paper/fig3/crop_organoid_umap.png', width=11, height=7, unit='cm')

ggplot(meta, aes(-UMAP1, -UMAP2, color=lineage)) +
  geom_point(size=0.5) +
  theme_void() + no_legend() +
  scale_color_manual(values=lineage_colors)
ggsave('plots/paper/fig3/crop_lineage_umap.png', width=11, height=7, unit='cm')

feature_plot(crop, features=c('EMX1', 'DLX5', 'LHX5', 'FOXG1'), order=T, pt.size=0.1) &
  scale_color_gradientn(colors=gyylgnbu())
ggsave('plots/paper/fig3/crop_feature_umap.png', width=18, height=15, unit='cm')

guide_dim_plot(crop, guide_assay = 'guide_assignments_cons')
ggsave('plots/paper/v2/fig3/guide_dim_plot_cons_all.png', width=20, height=12)


guide_meta <- crop_guide_counts %>% as_tibble(rownames='cell') %>%
    pivot_longer(!cell, names_to = 'guide', values_to = 'count') %>%
    group_by(cell) %>%
    arrange(desc(count)) %>%
    filter(row_number()==1) %>%
    right_join(as_tibble(crop_all@meta.data, rownames='cell')) %>%
    inner_join(as_tibble(crop_all[['umap']]@cell.embeddings, rownames='cell')) %>%
    mutate(target_gene=str_replace(guide, '(.+)_\\d', '\\1')) %>%
    filter(RSSUMAP_1<10) %>%
    arrange(!(count==0))

ggplot(arrange(guide_meta, !is.na(count)), aes(-RSSUMAP_1, -RSSUMAP_2, fill=target_gene)) +
    geom_point(size=3, shape=21, color='white', stroke=0.2) +
    theme_void() + no_legend() +
    scale_fill_manual(values=target_colors)
ggsave('~/projects/early/plots/paper/fig3/guide_confetti_p3_umap.png', width=16, height=12)




#### Enrichment in lineages ####
lineage_cmh_enrich <- test_guide_enrichment(
    crop,
    test_groups = 'lineage',
    min_counts_per_group = 10,
    groups = 'org_group',
    method = 'cmh'
)

lineage_fisher_enrich <- test_guide_enrichment(
    crop,
    test_groups = 'lineage',
    min_counts_per_group = 10,
    groups = 'org_group',
    method = 'fisher'
)


plot_df <- lineage_fisher_enrich %>%
    mutate(target_gene=str_replace(x, '(.+)-\\d', '\\1')) %>%
    mutate(padj=p.adjust(pval, method='fdr'))

ggplot(plot_df, aes(x, log_odds_ratio, color=group)) +
    geom_point(stat='identity') +
    geom_hline(yintercept = 0) +
    facet_grid(y~target_gene, scales='free')


#### Full plot ####

plot_df <- lineage_fisher_enrich %>%
    mutate(target_gene=str_replace(x, '(.+)-\\d', '\\1')) %>%
    filter(target_gene!='DUMMY')
ggplot(plot_df, aes(x, log_odds_ratio, color=group, size=-log10(pval), alpha=pval<0.05)) +
    geom_point(stat='identity') +
    geom_hline(yintercept = 0) +
    facet_grid(y~target_gene, scales='free') +
    no_x_text() + labs(x='guide') +
    ggtitle('Enrichment in lineages per organoid')
ggsave('plots/CROP/CROP_3/crop_3_all/lineage_enrichment_per_organoid.png', width=16, height=6)



#### Lollipop plot ####
### Here we filter guides that are significant in an individual organoid per fisher test
ftest_org_const_df <- lineage_fisher_enrich %>%
    dplyr::filter(pval<5e-2 & log_odds_ratio!=0) %>%
    dplyr::group_by(x, y) %>%
    dplyr::summarize(
        const_dir=max(sign(log_odds_ratio))==min(sign(log_odds_ratio)),
        mean_logodds_sig  = mean(log_odds_ratio)
    )

### Here we check if the enrichment is
# 1) consistent across organoids
# 2) shared by more than one guide targeting the same gene
ftest_sum_df <- lineage_fisher_enrich %>%
    # dplyr::filter(log_odds_ratio!=0) %>%
    dplyr::group_by(x, y) %>%
    dplyr::summarize(
        mean_logodds = mean(log_odds_ratio),
        sd_logodds=sd(log_odds_ratio),
    ) %>%
    ungroup() %>%
    dplyr::left_join(ftest_org_const_df) %>%
    dplyr::left_join(lineage_cmh_enrich) %>%
    dplyr::mutate(direction=sign(mean_logodds_sig), target_gene=str_replace(x, '(.+)-\\d', '\\1')) %>%
    dplyr::group_by(y, direction, target_gene) %>%
    dplyr::mutate(guide_const=all(const_dir) & n() > 1) %>%
    dplyr::mutate(show_label=guide_const&target_gene%in%c('GLI3', 'TBR1', 'ST18', 'SOX5', 'ZEB2', 'BACH2', 'FOXN4')) %>%
    mutate(padj=p.adjust(pval, method='fdr')) %>%
    mutate(
        signed_pclip=direction*-log10(ifelse(padj>1e-10, padj, 1e-10))) %>%
    ungroup() %>% group_by(y) %>%
    group_split() %>%
    map(function(x){
        x <- x %>% arrange(desc(signed_pclip)) %>%
            mutate(x=factor(x, levels=unique(.$x)))
    })


strokesize <- 0.2
size_range <- c(0.2, 2)
lolli_new <- function(df){
  df <- filter(df, !is.na(signed_pclip))
  ggplot(df, aes(x, signed_pclip, col=abs(signed_pclip)>4, fill=abs(signed_pclip)>4, size=abs(signed_pclip))) +
    geom_bar(aes(alpha=const_dir), stat='identity', color=NA, size=0.1, width=0.1) +
    geom_hline(yintercept=0, color='darkgray', size=strokesize) +
    geom_hline(yintercept=-4, color='darkgray', linetype='dotted', size=strokesize) +
    geom_hline(yintercept=4, color='darkgray', linetype='dotted', size=strokesize) +
    # Not consistent
    geom_point(data=df, shape=1, stroke=strokesize, alpha=0.8) +
    # Consistent across guides -> no alpha
    geom_point(data=dplyr::filter(df, const_dir & !guide_const), shape=21, stroke=0.1, alpha=0.3) +
    # Consistent across organoids -> with alpha
    geom_point(data=dplyr::filter(df, guide_const), shape=21, stroke=0.1, alpha=0.95) +
    # geom_point(data=dplyr::filter(ftest_sum_df, dir_const_sig), shape=1, stroke=strokesize, alpha=0.95) +
    scale_color_manual(values=c('darkgrey', 'black')) +
    scale_fill_manual(values=c('darkgrey', 'black')) +
    scale_alpha_manual(values=c(0.3,0.9)) +
    no_legend() + theme_rangeframe() + article_text() +
    scale_x_discrete(expand=c(0.05,0)) +
    # scale_y_continuous(expand=c(0.5,0.5), limits=c(-10,10), breaks=seq(-10,10,10)) +
    scale_axis_rangeframe() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(t=0,r=0,b=0,l=0)
    ) +
    # Labels up for enrichment
    geom_text(
      data=dplyr::filter(df, direction==1, guide_const),
      mapping = aes(label=x), y=13, col='black', size=5/ggplot2:::.pt,
      nudge_y = 30,
      nudge_x = 0,
      angle = 90,
      vjust = 0.5,
      hjust = 0
    ) +
    # Labels down for depletion
    geom_text(
      data=dplyr::filter(df, direction==-1, guide_const),
      mapping = aes(label=x), y=-13, col='black', size=5/ggplot2:::.pt,
      nudge_y = -5,
      nudge_x = 0,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ) +
    labs(y=expression(-log[10](q)))
}


plots <- map(ftest_sum_df, lolli_new)
wrap_plots(plots, ncol=1)

wrap_plots(plots[1:2], ncol=1) &
  scale_y_continuous(expand=c(0.5,0.5), limits=c(-10,10), breaks=seq(-10,10,10)) &
  scale_size_continuous(range=c(0.1,1.5))
ggsave('plots/paper/v2/fig3/guide_ge_ctx_bw_lolli.pdf', width=7.2, height=6, unit='cm')


#### Heatmap with consistent guides ####

consistent_guides <- ftest_sum_df %>%
  bind_rows() %>%
  filter(guide_const) %>%
  {.$x} %>% unique() %>% sort()

cons_guide_counts <- crop[['guide_assignments']]@data[consistent_guides, ]
gene_names <- stringr::str_replace(rownames(cons_guide_counts), '(.+)-\\d+', '\\1')
cons_target_counts <- cons_guide_counts %>%
  aggregate_matrix(groups=gene_names, fun='max')
colnames(cons_target_counts) <- colnames(cons_guide_counts)
crop[['target_genes_cons']] <- cons_target_counts %>%
  CreateAssayObject()
crop[['guide_assignments_cons']] <- cons_guide_counts %>%
  CreateAssayObject()

guide_dim_plot(crop, guide_assay='guide_assignments_cons')
ggsave('plots/paper/fig3/guide_dim_plot_umap.png', width=60, height=38, unit='cm')


consistent_cmh_enrich <- test_guide_enrichment(
  crop,
  guide_assay = 'target_genes_cons',
  test_groups = 'lineage',
  min_counts_per_group = 10,
  groups = 'org_group',
  method = 'cmh'
)

consistent_cmh_enrich$padj <- consistent_cmh_enrich$pval %>% p.adjust(method='bonferroni')


const_df <- ftest_sum_df %>%
  bind_rows() %>%
  dplyr::filter(guide_const, padj<1e-2) %>%
  dplyr::select(y, direction, 'x'=target_gene) %>%
  group_by(x,y) %>%
  dplyr::summarise(const=T, guide_support=n()) %>%
  filter(guide_support>1)

# Calculate how many guides were involved in the analysis
nguide_df <- ftest_sum_df %>%
  bind_rows() %>%
  dplyr::filter(guide_const) %>%
  dplyr::group_by(target_gene) %>%
  dplyr::summarise(n_guide=length(unique(x))) %>%
  dplyr::rename('x'=target_gene)

plot_df <- consistent_cmh_enrich %>%
  filter(y!='other') %>%
  left_join(const_df) %>%
  left_join(nguide_df) %>%
  group_by(x) %>% filter(any(!is.na(guide_support)))


bind_rows('cmh'=lineage_cmh_enrich, 'fisher'=lineage_fisher_enrich, 'cmh_consistent'=plot_df, .id='test_method') %>%
    write_tsv('data/CROP_3/crop_3_all/enrichment/all_enrichment_results.tsv')


plot_df <- read_tsv('data/CROP_3/crop_3_all/enrichment/all_enrichment_results.tsv') %>%
    filter(test_method=='cmh_consistent')


lodds_mat <- plot_df %>%
  dplyr::select(x, y, log_odds_ratio) %>%
  tidyr::spread(y, log_odds_ratio) %>%
  as_matrix()

issig_mat <- plot_df %>%
  dplyr::select(x, y, guide_support) %>%
  tidyr::spread(y, guide_support) %>%
  as_matrix() %>% {.[is.na(.)] <- 0; .} %>% binarize()

# gene_clust <- (lodds_mat * issig_mat) %>% dist() %>% hclust(method='ward.D2')
gene_clust <- lodds_mat %>% dist() %>% hclust(method='ward.D2')

tf_order <- c('GLI3', 'TBR1', 'ST18', 'BCL11B', 'PAX6', 'NEUROD1', 'SOX1', 'ZEB2', 'FOXN4', 'SOX9', 'NEUROD6', 'BACH2', 'SOX5', 'MYT1L', 'HES1', 'HOPX')
unique(plot_df$x) %in% tf_order


plot_df <- plot_df %>%
  mutate(x=factor(x, levels=rev(tf_order)),
         y=factor(y, levels=c('ctx', 'ge', 'nt')))


clim <- max(abs(plot_df$log_odds_ratio))


p1 <- ggplot(plot_df, aes(y, x, fill=log_odds_ratio)) +
  geom_tile() +
  # geom_point(mapping=aes(shape=factor(guide_support)), size=0.5, alpha=0.99, color='#263238') +
  geom_point(data=filter(plot_df, !is.na(guide_support)), size=1.5, alpha=0.99, color='#263238', shape=21, fill='darkgrey', stroke=0.1) +
  scale_fill_gradientn(colors=rev(pals::brewer.rdylbu(100)), limits=c(-1.5, 1.5)) +
  scale_x_discrete(expand=c(0,0), labels=c('Dorsal', 'Ventral', 'N.t.')) +
  scale_y_discrete(expand=c(0,0)) +
  no_label() +
  no_legend() +
  article_text() +
  theme(
    plot.margin = unit(rep(0, 4), 'lines'),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size=0.3)
  )
p1

p2 <- ggplot(plot_df, aes(1, x, fill=factor(n_guide))) +
  geom_tile(color='black', size=0.1) +
  no_label() +
  no_text() +
  scale_fill_manual(values=c('gray', 'black')) +
  no_legend() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(
    plot.margin = unit(rep(0, 4), 'lines'),
    panel.border = element_blank()
  )

p1 + p2 + plot_layout(widths=c(5,0.5))

ggsave('plots/paper/v2/fig3/guide_consistent_heatmap.pdf', width=3, height=4, unit='cm')


print_scale(pals::brewer.rdylbu(100))
ggsave('plots/paper/fig3/brewer_rdylbu.pdf')



#### DE globally ####
crop_neurons <- subset(crop, celltype%in%c('lge_in', 'ctx_ex', 'mesen_in'))
crop_wt_mat <- t(GetAssayData(crop_neurons, assay = 'RNA', slot = 'data'))
ko_prob <- t(crop_neurons[['ko_prob_global']]@data)
cells_use <- intersect(rownames(crop_wt_mat), rownames(ko_prob))

covariate_names <- c('nFeature_RNA', 'organoid', 'celltype')
C <- crop@meta.data[cells_use, covariate_names]

target_genes <- str_replace(colnames(ko_prob), '(.+)-\\d', '\\1')
X_new <- t(aggregate_matrix(t(as.matrix(ko_prob)), groups = target_genes, fun=colMaxs))
rownames(X_new) <- rownames(ko_prob)
XC_new <- cbind(C, X_new)

model_formula <- reformulate(termlabels = c(covariate_names, colnames(X_new)))
XC_model_new <- model.matrix(model_formula, XC_new)
test_terms <- which(colnames(XC_model_new) %in% colnames(X_new))

prc_expr <- colMeans(crop_wt_mat[cells_use, ]>0)
expr_genes <- prc_expr[prc_expr>0.05 & !str_detect(names(prc_expr), '^(MT-|RP|MT\\.)') & names(prc_expr)%in%good_genes]

Y <- crop_wt_mat[cells_use, names(expr_genes)]

refit_de <- glm_de(XC_model_new, Y, term = test_terms)
refit_de$padj <- p.adjust(refit_de$pval, method = 'fdr')

refit_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_neurons_de.tsv')
# refit_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_neurons_de.tsv')





#### Number of DEG vs composition change ####

n_deg <- refit_de %>%
  filter(padj<1e-4) %>%
  group_by(name) %>%
  summarize(ndeg=n(), target_gene=name) %>%
  distinct()

n_cells <- colSums(X_new>0) %>%
  enframe('target_gene', 'ncells')

enrich_vs_de <- consistent_cmh_enrich %>%
  group_by(x) %>%
  filter(log_odds_ratio==max(log_odds_ratio)) %>%
  mutate(enrich_score=abs(log_odds_ratio)) %>%
  mutate(target_gene=x) %>%
  right_join(n_deg) %>%
  right_join(n_cells) %>%
  mutate(enrich_score=ifelse(is.na(enrich_score), 0, enrich_score))


ggplot(enrich_vs_de, aes(ndeg, enrich_score, label=target_gene, size=log10(ncells))) +
  geom_point(shape=21, fill='black', color='white', stroke=0.2) +
  geom_text_repel(size=5/ggplot2::.pt, max.overlaps = 100) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base=10), breaks=c(0,10,100,1000), limits=c(0,1000)) +
  scale_y_continuous(limits=c(0,1.5)) +
  scale_size_continuous(range=c(0.5,3)) +
  theme_rangeframe() + scale_axis_rangeframe() +
  no_legend() + article_text() +
  labs(x='# DEG (state change)', y='Composition change')

ggsave('plots/paper/v2/fig3/deg_vs_comp_change_scatter.pdf', width=4.6, height=4.5, units='cm')



print_scale(pals::brewer.rdbu(100))
ggsave('plots/paper/v2/fig3/rdbu.pdf', width=4.6, height=4.5, units='cm')


#### Get DE in neurons ####
## CROP ####
#### Ventral ####
crop_ge <- subset(crop, celltype=='lge_in')

crop_wt_mat <- t(GetAssayData(crop_ge, assay = 'RNA', slot = 'data'))
cells_use <- intersect(rownames(crop_wt_mat), rownames(crop_guide_counts))

X <- binarize(crop_guide_counts[cells_use, ])
target_genes <- str_replace(colnames(X), '(.+)_\\d', '\\1')

covariate_names <- c('nFeature_RNA', 'organoid')
C <- crop_ge@meta.data[cells_use, covariate_names]
model_formula <- reformulate(termlabels = c(covariate_names, colnames(X)))
XC <- cbind(C, X)
XC_model <- model.matrix(model_formula, XC)

var_feats <- variable_features(crop_wt_mat[cells_use, ])
var_genes <- var_feats %>%
  dplyr::filter(!str_detect(gene, '^(MT-|RP|MT\\.)')) %>%
  dplyr::top_n(500, vst_variance_standardized) %>%
  {.$gene}
Y <- crop_wt_mat[cells_use, var_genes]

ko_probs <- infer_ko_probabilities(X, Y, C, alpha=0.5, target_genes=target_genes, parallel=T)

target_genes <- str_replace(colnames(ko_probs$ko_probabilities), '(.+)_\\d', '\\1')
X_new <- t(aggregate_matrix(t(as.matrix(ko_probs$ko_probabilities)), groups = target_genes, fun=colMaxs))
rownames(X_new) <- rownames(ko_probs$ko_probabilities)

# X_new <- t(aggregate_matrix(t(X), groups = target_genes, fun=colMaxs))
# rownames(X_new) <- rownames(X)

XC_new <- cbind(C, X_new)
model_formula <- reformulate(termlabels = c(covariate_names, colnames(X_new)))
XC_model_new <- model.matrix(model_formula, XC_new)
test_terms <- which(colnames(XC_model_new) %in% colnames(X_new))

prc_expr <- colMeans(crop_wt_mat[cells_use, ]>0)
expr_genes <- prc_expr[prc_expr>0.05 & !str_detect(names(prc_expr), '^(MT-|RP|MT\\.)') & names(prc_expr)%in%good_genes]

Y <- crop_wt_mat[cells_use, names(expr_genes)]

refit_de <- glm_de(XC_model_new, Y, term = test_terms)

perturb_de <- refit_de %>%
  filter(name=='GLI3')
perturb_de$padj <- p.adjust(perturb_de$pval, method = 'fdr')
refit_de$padj <- p.adjust(refit_de$pval, method = 'fdr')

# perturb_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_gli3_ge_in_de.tsv')
# refit_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ge_in_de.tsv')
refit_de_new <- refit_de

refit_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ge_in_de.tsv')


#### Ventral NPCs ####
crop_ge <- subset(crop, celltype=='ge_npc')

crop_wt_mat <- t(GetAssayData(crop_ge, assay = 'RNA', slot = 'data'))
cells_use <- intersect(rownames(crop_wt_mat), rownames(crop_guide_counts))

X <- binarize(crop_guide_counts[cells_use, ])
target_genes <- str_replace(colnames(X), '(.+)_\\d', '\\1')

covariate_names <- c('nFeature_RNA', 'organoid', 'Phase')
C <- crop_ge@meta.data[cells_use, covariate_names]
model_formula <- reformulate(termlabels = c(covariate_names, colnames(X)))
XC <- cbind(C, X)
XC_model <- model.matrix(model_formula, XC)

var_feats <- variable_features(crop_wt_mat[cells_use, ])
var_genes <- var_feats %>%
  dplyr::filter(!str_detect(gene, '^(MT-|RP|MT\\.)')) %>%
  dplyr::top_n(500, vst_variance_standardized) %>%
  {.$gene}
Y <- crop_wt_mat[cells_use, var_genes]

ko_probs <- infer_ko_probabilities(X, Y, C, alpha=0.5, target_genes=target_genes, parallel=T)

target_genes <- str_replace(colnames(ko_probs$ko_probabilities), '(.+)_\\d', '\\1')
X_new <- t(aggregate_matrix(t(as.matrix(ko_probs$ko_probabilities)), groups = target_genes, fun=colMaxs))
rownames(X_new) <- rownames(ko_probs$ko_probabilities)

# X_new <- t(aggregate_matrix(t(X), groups = target_genes, fun=colMaxs))
# rownames(X_new) <- rownames(X)

XC_new <- cbind(C, X_new)
model_formula <- reformulate(termlabels = c(covariate_names, colnames(X_new)))
XC_model_new <- model.matrix(model_formula, XC_new)
test_terms <- which(colnames(XC_model_new) %in% colnames(X_new))

prc_expr <- colMeans(crop_wt_mat[cells_use, ]>0)
expr_genes <- prc_expr[prc_expr>0.05 & !str_detect(names(prc_expr), '^(MT-|RP|MT\\.)') & names(prc_expr)%in%good_genes]

Y <- crop_wt_mat[cells_use, names(expr_genes)]

refit_de <- glm_de(XC_model_new, Y, term = test_terms)

perturb_de <- refit_de %>%
  filter(name=='GLI3')
perturb_de$padj <- p.adjust(perturb_de$pval, method = 'fdr')
refit_de$padj <- p.adjust(refit_de$pval, method = 'fdr')

perturb_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_gli3_ge_npc_de.tsv')
refit_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ge_npc_de.tsv')

refit_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ge_npc_de.tsv')

ggplot(filter(refit_de, padj<0.05), aes(coef, -log10(padj), label=gene)) +
  geom_text(size=2) +
  facet_wrap(~name, scales='free')


#### Dorsal ####
crop_ctx <- subset(crop, celltype=='ctx_ex')

crop_wt_mat <- t(GetAssayData(crop_ctx, assay = 'RNA', slot = 'data'))
cells_use <- intersect(rownames(crop_wt_mat), rownames(crop_guide_counts))

X <- binarize(crop_guide_counts[cells_use, ])
target_genes <- str_replace(colnames(X), '(.+)_\\d', '\\1')

covariate_names <- c('nFeature_RNA', 'organoid')
C <- crop_ctx@meta.data[cells_use, covariate_names]
model_formula <- reformulate(termlabels = c(covariate_names, colnames(X)))
XC <- cbind(C, X)
XC_model <- model.matrix(model_formula, XC)

var_feats <- variable_features(crop_wt_mat[cells_use, ])
var_genes <- var_feats %>%
  dplyr::filter(!str_detect(gene, '^(MT-|RP|MT\\.)')) %>%
  dplyr::top_n(500, vst_variance_standardized) %>%
  {.$gene}
Y <- crop_wt_mat[cells_use, var_genes]

ko_probs <- infer_ko_probabilities(X, Y, C, alpha=0.5, target_genes=target_genes, parallel=T)

target_genes <- str_replace(colnames(ko_probs$ko_probabilities), '(.+)_\\d', '\\1')
X_new <- t(aggregate_matrix(t(as.matrix(ko_probs$ko_probabilities)), groups = target_genes, fun=colMaxs))
rownames(X_new) <- rownames(ko_probs$ko_probabilities)

XC_new <- cbind(C, X_new)
model_formula <- reformulate(termlabels = c(covariate_names, colnames(X_new)))
XC_model_new <- model.matrix(model_formula, XC_new)
test_terms <- which(colnames(XC_model_new) %in% colnames(X_new))

prc_expr <- colMeans(crop_wt_mat[cells_use, ]>0)
expr_genes <- prc_expr[prc_expr>0.05 & !str_detect(names(prc_expr), '^(MT-|RP|MT\\.)') & names(prc_expr)%in%good_genes]
Y <- crop_wt_mat[cells_use, names(expr_genes)]

refit_de <- glm_de(XC_model_new, Y, term = test_terms)
perturb_de <- refit_de %>%
  filter(name=='GLI3')
perturb_de$padj <- p.adjust(perturb_de$pval, method = 'fdr')
refit_de$padj <- p.adjust(refit_de$pval, method = 'fdr')

perturb_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_gli3_ctx_ex_de.tsv')
refit_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ctx_ex_de.tsv')

refit_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ctx_ex_de.tsv')

ggplot(filter(refit_de, padj<0.05), aes(coef, -log10(padj), label=gene)) +
  geom_text(size=2) +
  facet_wrap(~name, scales='free')


#### Dorsal NPC ####
crop_ctx_npc <- subset(crop, celltype=='ctx_npc')

crop_wt_mat <- t(GetAssayData(crop_ctx_npc, assay = 'RNA', slot = 'data'))
cells_use <- intersect(rownames(crop_wt_mat), rownames(crop_guide_counts))

X <- binarize(crop_guide_counts[cells_use, ])
target_genes <- str_replace(colnames(X), '(.+)_\\d', '\\1')

covariate_names <- c('nFeature_RNA', 'organoid', 'Phase')
C <- crop_ctx_npc@meta.data[cells_use, covariate_names]
model_formula <- reformulate(termlabels = c(covariate_names, colnames(X)))
XC <- cbind(C, X)
XC_model <- model.matrix(model_formula, XC)

var_feats <- variable_features(crop_wt_mat[cells_use, ])
var_genes <- var_feats %>%
  dplyr::filter(!str_detect(gene, '^(MT-|RP|MT\\.)')) %>%
  dplyr::top_n(500, vst_variance_standardized) %>%
  {.$gene}
Y <- crop_wt_mat[cells_use, var_genes]

ko_probs <- infer_ko_probabilities(X, Y, C, alpha=0.5, target_genes=target_genes, parallel=T)

target_genes <- str_replace(colnames(ko_probs$ko_probabilities), '(.+)_\\d', '\\1')
X_new <- t(aggregate_matrix(t(as.matrix(ko_probs$ko_probabilities)), groups = target_genes, fun=colMaxs))
rownames(X_new) <- rownames(ko_probs$ko_probabilities)

XC_new <- cbind(C, X_new)
model_formula <- reformulate(termlabels = c(covariate_names, colnames(X_new)))
XC_model_new <- model.matrix(model_formula, XC_new)
test_terms <- which(colnames(XC_model_new) %in% colnames(X_new))

prc_expr <- colMeans(crop_wt_mat[cells_use, ]>0)
expr_genes <- prc_expr[prc_expr>0.05 & !str_detect(names(prc_expr), '^(MT-|RP|MT\\.)') & names(prc_expr)%in%good_genes]
Y <- crop_wt_mat[cells_use, names(expr_genes)]

refit_de <- glm_de(XC_model_new, Y, term = test_terms)
perturb_de <- refit_de %>%
  filter(name=='GLI3')
perturb_de$padj <- p.adjust(perturb_de$pval, method = 'fdr')
refit_de$padj <- p.adjust(refit_de$pval, method = 'fdr')

perturb_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_gli3_ctx_npc_de.tsv')
refit_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ctx_npc_de.tsv')

refit_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ctx_npc_de.tsv')

ggplot(filter(refit_de, padj<0.05), aes(coef, -log10(padj), label=gene)) +
  geom_text(size=2) +
  facet_wrap(~name, scales='free')


### Compare ctx and GE ####
ge_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ge_in_de.tsv')
ctx_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ctx_ex_de.tsv')

compare_de <- ge_de %>%
  inner_join(ctx_de, by=c('gene', 'name')) %>%
  filter(padj.x<0.05 | padj.y<0.05)

p1 <- ggplot(compare_de, aes(coef.x, coef.y, label=gene)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_text(size=2) +
  facet_wrap(~name)
p1

p1 <- ggplot(compare_de, aes(coef.x, coef.y, label=gene)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_pointdensity(size=0.2, adjust=0.5) +
  geom_smooth(method='lm') +
  facet_wrap(~name, scales='free')
p1



### Combine results ####
ge_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ge_in_de.tsv')
ctx_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ctx_ex_de.tsv')
ge_npc_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ge_npc_de.tsv')
ctx_npc_de <- read_tsv('data/CROP_3/crop_3_all/diff_expression/crop_ctx_npc_de.tsv')

crop_de <- bind_rows('ventral_neurons'=ge_de, 'dorsal_neurons'=ctx_de, 'ventral_npcs'=ge_npc_de, 'dorsal_npcs'=ctx_npc_de) %>%
  filter(padj<0.05)

crop_de %>% write_tsv('data/CROP_3/crop_3_all/diff_expression/crop_all_de.tsv')



#### Sub GRNs for selected target genes ####

npc_idx <- rnatac$pseudotime_ranks>0.1 & rnatac$pseudotime_ranks<0.55 & rnatac$lineage%in%c('telencephalon', 'ge', 'ctx')
npc_expr_frac <- rowMeans(GetAssayData(rnatac, assay='RNA')[, npc_idx]>0)
npc_genes <- names(npc_expr_frac)[npc_expr_frac>0.05]
npc_targets <- intersect(npc_genes, target_names)


module_peaks <- modules %>%
    separate_rows(peaks, sep=';')

grn_regions <- module_peaks %>% pull(peaks) %>% str_replace_all('_', '-')
grn_peaks <- unique(grnreg2peaks[grn_regions])
grn_peaks_acc <- rnatac@assays$peaks@counts[grn_peaks, ]


rnatac$this <- rnatac$pseudotime_ranks>0.1 & rnatac$pseudotime_ranks<0.55
dim_plot(rnatac, group.by='this')

branch_clusters <- cluster_meta %>%
    filter(pseudotime_ranks>0.3 & pseudotime_ranks<0.9 & lineage%in%c('telencephalon', 'ge', 'ctx')) %>% pull(name)

ggplot(cluster_meta, aes(FR1, FR2, color=name%in%branch_clusters, label=name)) +
    geom_text()

highres_acc <- aggregate_matrix(t(grn_peaks_acc > 0), groups = rnatac$highres_clusters)
is_acc <- highres_acc[branch_clusters, ] > 0.1
acc_peaks <- colnames(is_acc)[colMeans2(is_acc)>0]
acc_ranges <- acc_peaks %>% StringToGRanges()
acc_regions <- names(grnreg2peaks)[match(acc_peaks, grnreg2peaks)] %>% str_replace_all('-', '_')

branch_modules <- module_peaks %>%
    filter(peaks%in%acc_regions) %>%
    group_by(tf, target, estimate, corr) %>% summarize(peaks=paste(peaks, collapse=';'))

# branch_modules %>% write_tsv('~/projects/early/data/RNA_ATAC/grn/networks/v5_branchpoint_modules.tsv')

#### Correlation to lineage ####
hc_meta_use <- cluster_meta %>%
    filter(pseudotime_ranks>0.1)

match_idx <- match(hc_meta_use$name, rownames(rnatac@misc$cluster_summaries$RNA))
ctx_rna_corr <- cor(hc_meta_use$to_ctx_ranks, as.matrix(rnatac@misc$cluster_summaries$RNA[match_idx, ]), method = 'pearson')
ge_rna_corr <- cor(hc_meta_use$to_ge_ranks, as.matrix(rnatac@misc$cluster_summaries$RNA[match_idx, ]), method = 'pearson')
nt_rna_corr <- cor(hc_meta_use$to_nt_ranks, as.matrix(rnatac@misc$cluster_summaries$RNA[match_idx, ]), method = 'pearson')


lineage_probs <- do.call(rbind, list(ctx_rna_corr, ge_rna_corr, nt_rna_corr)) %>% t() %>%
    as_tibble(rownames='gene') %>% {colnames(.)[-1] <- c('ctx_corr', 'ge_corr', 'nt_corr'); .} %>%
    filter(!is.na(ge_corr))


gene_meta_ <- gene_meta %>%
    select(-c(ctx_corr, ge_corr, nt_corr)) %>%
    left_join(lineage_probs)


module_meta <- branch_modules %>%
    filter(target%in%npc_genes) %>%
    inner_join(gene_meta_, by=c('target'='gene')) %>%
    mutate(dir=sign(estimate))


module_lineage_hr <- module_meta %>%
    group_by(tf, dir) %>%
    summarize(
        ctx_sum = sum(ctx_corr),
        ge_sum = sum(ge_corr),
        nt_sum = sum(nt_corr),
        ctx_mean = mean(ctx_corr),
        ge_mean = mean(ge_corr),
        nt_mean = mean(nt_corr)
    ) %>%
    group_by(tf) %>%
    mutate(
        ctx_diff = ctx_mean[2] - ctx_mean[1],
        ge_diff = ge_mean[2] - ge_mean[1],
        nt_diff = nt_mean[2] - nt_mean[1]
    )



g1_module <- filter(module_meta, tf=='GLI3')

ggplot(g1_module, aes(factor(dir), ctx_corr)) +
    geom_boxplot()

ggplot(g1_module, aes(factor(dir), ge_corr)) +
    geom_boxplot()

ggplot(g1_module, aes(factor(dir), nt_corr)) +
    geom_boxplot()


g2_module <- filter(module_meta, tf=='HES1')


ggplot(g2_module, aes(factor(dir), ctx_corr)) +
    geom_boxplot()

ggplot(g2_module, aes(factor(dir), ge_corr)) +
    geom_boxplot()

ggplot(g2_module, aes(factor(dir), nt_corr)) +
    geom_boxplot()


target_intersect <- intersect(g2_module$target, g1_module$target)
target_union <- union(g2_module$target, g1_module$target)

inter_modules <- filter(module_meta, target%in%target_union, tf%in%c('HES1', 'GLI3')) %>%
    mutate(
        group=case_when(
            target%in%target_intersect ~ 'intersect',
            target%in%g2_module$target ~ 'HES1',
            target%in%g1_module$target ~ 'GLI3'
        )
    )

ggplot(inter_modules, aes(factor(dir), ge_corr, fill=group)) +
    facet_grid(~tf) +
    geom_boxplot()

ggplot(inter_modules, aes(factor(dir), ctx_corr, fill=group)) +
    facet_grid(~tf) +
    geom_boxplot()

ggplot(inter_modules, aes(factor(dir), nt_corr, fill=group)) +
    facet_grid(~tf) +
    geom_boxplot()


#### Plot subGRN ####
gli3_module <- filter(module_meta, tf=='GLI3')
gli3_graph <- gli3_module %>%
    as_tbl_graph() %N>%
    left_join(gli3_module, by=c('name'='target'))


ggraph(gli3_graph, x=ctx_corr, y=as.numeric(name=='GLI3')) +
    # geom_edge_diagonal(mapping=aes(color=sign(estimate))) +
    geom_edge_diagonal(mapping=aes(color=estimate)) +
    geom_node_point() +
    # geom_node_text(aes(label=name), size=8/ggplot2::.pt, angle=90) +
    scale_edge_color_gradientn(colors=rev(pals::ocean.curl(100)[-c(42:58)]), trans='pseudo_log', limits=c(-1.3, 1.3)) +
    # scale_edge_color_gradientn(colors=c('#B94D60', '#70AD8C'), limits=c(-1.1,1.1)) +
    scale_x_continuous(na.value = 0) +
    # no_legend() +
    theme_void()
ggsave('plots/paper/v2/fig3/gli3_subnet.pdf')


plot_df <- gli3_module
p1 <- ggplot(plot_df, aes(ctx_corr, factor(dir), label=target)) +
    geom_boxplot(fill=lineage_colors['ctx'], alpha=0.5, size=0.2, outlier.size = 0.2) +
    # geom_point(aes(alpha=target=='ARID1B')) +
    scale_x_continuous(limits=c(-0.75, 0.75)) +
    theme_void()

p2 <- ggplot(plot_df, aes(ge_corr, factor(dir), label=target)) +
    geom_boxplot(fill=lineage_colors['ge'], alpha=0.5, size=0.2, outlier.size = 0.2) +
    # geom_text() +
    scale_x_continuous(limits=c(-0.75, 0.75), breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_y_text() +
    theme(
        axis.line.y = element_blank()
    ) +
    no_label()

p1 / p2 & no_margin()
ggsave('plots/paper/v2/fig3/gli3_subnet_box.pdf', width=5.2, height=1.5, unit='cm')




hes1_module <- filter(module_meta, tf=='HES1')
hes1_graph <- hes1_module %>%
    as_tbl_graph() %N>%
    left_join(hes1_module, by=c('name'='target'))


ggraph(hes1_graph, x=ctx_corr, y=as.numeric(name=='HES1')) +
    # geom_edge_diagonal(mapping=aes(color=sign(estimate))) +
    geom_edge_diagonal(mapping=aes(color=estimate)) +
    geom_node_point() +
    # geom_node_text(aes(label=name), size=8/ggplot2::.pt, angle=90) +
    scale_edge_color_gradientn(colors=rev(pals::ocean.curl(100)[-c(42:58)]), trans='pseudo_log', limits=c(-1.3, 1.3)) +
    # scale_edge_color_gradientn(colors=c('#B94D60', '#70AD8C'), limits=c(-1.1,1.1)) +
    scale_x_continuous(na.value = 0) +
    no_legend() +
    theme_void()
ggsave('plots/paper/v2/fig3/hes1_subnet.pdf')


plot_df <- hes1_module
p1 <- ggplot(plot_df, aes(ctx_corr, factor(dir), label=target)) +
    geom_boxplot(fill=lineage_colors['ctx'], alpha=0.5, size=0.2, outlier.size = 0.2) +
    # geom_point(aes(alpha=target=='ARID1B')) +
    scale_x_continuous(limits=c(-0.4, 0.4)) +
    theme_void()

p2 <- ggplot(plot_df, aes(ge_corr, factor(dir), label=target)) +
    geom_boxplot(fill=lineage_colors['ge'], alpha=0.5, size=0.2, outlier.size = 0.2) +
    # geom_text() +
    scale_x_continuous(limits=c(-0.4, 0.4), breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_y_text() +
    theme(
        axis.line.y = element_blank()
    ) +
    no_label()

p1 / p2 & no_margin()
ggsave('plots/paper/v2/fig3/hes1_subnet_box.pdf', width=4, height=1.5, unit='cm')







