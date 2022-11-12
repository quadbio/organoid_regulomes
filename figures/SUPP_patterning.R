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
# rnatac %>% write_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2.1_srt.rds')


#### Read objects ####
nect <- read_rds('data/RNA_ATAC/subsets/RNA_nect_srt.rds')
nepi <- read_rds('data/RNA_ATAC/subsets/RNA_nepi_srt.rds')
branch <- read_rds('data/RNA_ATAC/subsets/RNA_branch_srt.rds')


#### Feature plots ####
goi <- c('HES1', 'HES2', 'HES3', 'HES4', 'HES5', 'HES6', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'DLL1', 'DLL3', 'DLL4', 'JAGN1', 'JAG2', 'JAG1', 'DLK1')
p1 <- feature_plot(nect, features=goi, ncol=length(goi), order=T)
p2 <- feature_plot(nepi, features=goi, ncol=length(goi), order=T)
p3 <- feature_plot(branch, features=goi, ncol=length(goi), order=T)
p4 <- feature_plot(rnatac, features=goi, ncol=length(goi), order=T)
p5 <- feature_plot(gli3_muo_use, features=goi, ncol=length(goi), order=T)


p1 / p2 / p3 / p4 / p5
ggsave('plots/RNA_ATAC/annotation/hes_notch_umap.png', width=40, height=12)

feature_plot(gli3_muo_use, features=goi, ncol=length(goi), order=T, split.by='is_KO')
ggsave('plots/RNA_ATAC/annotation/hes_notch_gli3_ko_umap.png', width=8, height=40)


## Select allowed genes 
mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl')
all_coding_genes <- getBM(attributes = c( 'hgnc_symbol'), filters = c('biotype'), values = list(biotype='protein_coding'), mart = mart)
good_genes <- all_coding_genes$hgnc_symbol
good_genes <- good_genes[!str_detect(good_genes, '^MT-|^RP|^HIST')]
cc_genes_all <- cc.genes.updated.2019 %>% 
    reduce(union) %>% 
    union(cc_genes)

feature_sets <- read_rds('data/gene_sets/RNA_feature_sets.rds')

gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')
gene_annot <- gene_annot[gene_annot$gene_name %in% feature_sets$grn_features]
peak_ranges <- StringToGRanges(rownames(rnatac@assays$peaks))

lineage_graph <- read_rds('data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')
lineage_tree <- read_rds('data/RNA_ATAC/lineages/cellrank_lineage_tree.rds')
hc_meta <- read_tsv('data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv') %>% 
    mutate(cluster=as.character(cluster))
cluster_meta <- as_tibble(lineage_graph)


umap_meta <- Reductions(rnatac, slot='umap')@cell.embeddings %>% 
    as_tibble(rownames='cell')

meta <- rnatac@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(umap_meta)

gene_meta <- read_tsv('data/gene_sets/gene_scores.tsv')

organizer_de <- read_tsv('/links/groups/treutlein/USERS/jfleck/projects/early/data/gene_sets/lamanno_organizer_de.tsv')
organizer_markers <- read_tsv('/links/groups/treutlein/USERS/jfleck/projects/early/data/gene_sets/lamanno_organizer_markers.txt', col_names=F)$X1
organizers_npcs <- read_rds('/links/groups/treutlein/USERS/jfleck/data/mousebrain_la_manno/lamanno_mouse_dev_brain_organizers.rds')

organizer_mean <- t(aggregate_matrix(t(organizers_npcs@assays$RNA@data), organizers_npcs$organizer))

dim_plot(organizers_npcs, group.by='organizer', label=T)


#### Embed early stage (neuroectoderm) ####

rnatac$nect <- rnatac$age>=7 & rnatac$age<12 & rnatac$pseudotime_ranks<0.15 & rnatac$pseudotime_ranks>0.05 & rnatac$lineage!='other'

dim_plot(rnatac, group.by='nect')
dim_plot(rnatac, group.by='lineage')
feature_plot(rnatac, features='pseudotime_ranks')

nect <- subset(rnatac, nect==TRUE)

nect@active.assay <- 'RNA'
nect <- FindVariableFeatures(nect, nfeatures=200)
var_feats <- nect@assays$RNA@meta.features %>% 
    as_tibble(rownames='gene') %>%
    filter(gene%in%good_genes) %>% 
    top_n(100, vst.variance.standardized) %>% 
    pull(gene) %>% unique()
VariableFeatures(nect) <- var_feats %>% union(organizer_markers) %>% setdiff(cc_genes_all)
# VariableFeatures(nect) <- organizer_markers %>% setdiff(cc_genes_all)
nect <- CellCycleScoring(nect, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
nect <- ScaleData(nect, vars.to.regress = c('S.Score', 'G2M.Score'))
nect <- nect %>% RunPCA() 

### -> no integration here

nect <- RunUMAP(
    object = nect,
    spread = 0.5,
    min.dist = 0.2,
    reduction = 'pca', dims = 1:10,
    reduction.name = 'umap',
    reduction.key = 'UMAP_'
)

dim_plot(nect, group.by=c('age', 'line')) 

feature_plot(nect, reduction='umap', features=c('KRT8', 'DLX5', 'GLI3'), order=T)


nect <- nect %>% FindNeighbors(reduction = 'pca', dims = 1:10) 
nect <- nect %>% FindClusters(resolution=1)
nect_cluster_expr <- nect@assays$RNA@data
nect_cluster_mean <- t(aggregate_matrix(t(nect_cluster_expr), nect$seurat_clusters))
nect_meta <- nect@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(Reductions(nect, 'umap')@cell.embeddings, rownames='cell'))

dim_plot(nect, group.by=c('age', 'line', 'lineage'))


#### Write results ####
nect %>% write_rds('data/RNA_ATAC/subsets/RNA_nect.srt')
nect_rna <- DietSeurat(nect, assays='RNA', dimreducs = c('umap', 'fa2', 'css', 'pca'))
nect_rna %>% SaveH5Seurat(filename = 'data/RNA_ATAC/subsets/RNA_nect.h5Seurat', overwrite=T)
Convert(source = 'data/RNA_ATAC/subsets/RNA_nect.h5Seurat', dest='h5ad', overwrite=T)



#### Correlate with organizers ####
genes_use <- intersect(organizer_markers, rownames(nect_cluster_mean))
organizer_cor <- cor(as.matrix(nect_cluster_mean[genes_use, ]), as.matrix(organizer_mean[genes_use, ]))
pheatmap::pheatmap(organizer_cor)

organizer_cor_df <- organizer_cor %>% 
    as_tibble(rownames='seurat_clusters') %>% pivot_longer(!seurat_clusters, names_to='organizer', values_to='corr') %>% 
    inner_join(nect_meta)



#### DE between clusters ####

nect_cluster_de <- de(nect, groups = 'seurat_clusters')

p1 <- dim_plot(nect, label=T)
p2 <- dim_plot(nect, group.by=c('line')) +
    scale_color_manual(values=line_colors)
p3 <- dim_plot(nect, group.by=c('age')) +
    scale_color_manual(values=age_colors)
p4 <- ggplot(organizer_cor_df, aes(UMAP_1, UMAP_2, color=corr)) +
    geom_point(size=0.3) +
    scale_color_gradientn(colors=gyrdpu()) +
    facet_wrap(~organizer) + theme_void()

p5 <- vulcano_plot(nect_cluster_de)

(p1 + p2 + p3) / p4 / p5 + plot_layout(heights=c(1,3,5))
ggsave('plots/RNA_ATAC/subsets/nect_overview.png', width=14, height=16)


feature_plot(nect, features=c('GATA2', 'GATA3', 'TFAP2A', 'TFAP2C', 'PAX2', 'DLX5', 'PAX3', 'PAX8', 'ST18'), order=T)
hcr_probes <- c('FOXG1', 'TCF7L2', 'ASCL1', 'WNT8B', 'FGF8', 'EMX1', 'DLX2', 'NKX2-1', 'TRH', 'SFRP2', 'NRG1', 
  'BARHL2','OTP', 'SHH', 'TBR1', 'MSX1', 'LMX1A', 'DCX','DLX5','HOXB2', 'SIX3','FGF19','IRX3',
  'ZIC1','WNT9A','BMP7','TAL2','RSPO2','RSPO3', 'FEZF1', 'SIX6', 'OTX2', 'FOXA2', 'SOX2', 'GBX2', 
  'GLI3', 'BMP7', 'CTIP2', 'BCL11B')

feature_plot(nect, features=hcr_probes, order=T)
ggsave('~/projects/early/plots/RNA_ATAC/subsets/nect_hcr.png', width=10, height=20)

feature_plot(et_nepi, features=hcr_probes, order=T)
ggsave('~/projects/early/plots/RNA_ATAC/subsets/nepi_hcr.png', width=10, height=20)

feature_plot(et_branch, features=hcr_probes, order=T)
ggsave('~/projects/early/plots/RNA_ATAC/subsets/branch_hcr.png', width=10, height=20)


#### Plots for figure ####

rnatac_meta <- rnatac@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(Reductions(rnatac, 'umap')@cell.embeddings, rownames='cell'))

plot_df <- rnatac_meta
ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=nect, alpha=nect)) +
    geom_point(size=2, shape=21, stroke=0.1) +
    scale_fill_manual(values=c('darkgrey', '#212121')) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/full_nect_umap.png', width=20, height=12, unit='cm')


plot_df <- nect_meta
ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=line)) +
    geom_point(size=3, shape=21, stroke=0.1) +
    scale_fill_manual(values=line_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nect_line_umap.png', width=15, height=10, unit='cm')

ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=factor(age))) +
    geom_point(size=3, shape=21, stroke=0.1) +
    scale_fill_manual(values=age_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nect_age_umap.png', width=15, height=10, unit='cm')

ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=pseudotime_ranks)) +
    geom_point(size=3, shape=21, stroke=0.1) +
    scale_fill_viridis(option='magma', limits=c(0,1), direction=-1) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nect_pt_umap.png', width=15, height=10, unit='cm')

cluster_umap_means <- Reductions(nect, 'umap')@cell.embeddings %>% 
    aggregate_matrix(nect$seurat_clusters) %>% as_tibble(rownames='seurat_clusters')
plot_df <- nect_meta 
ggplot(plot_df, aes(UMAP_1, UMAP_2, color=seurat_clusters)) +
    geom_point(size=1.5, alpha=0.5) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nect_clusters_umap.png', width=2.85*3, height=1.9*3, unit='cm')


ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=seurat_clusters)) +
    geom_label(data=cluster_umap_means, mapping=aes(label=seurat_clusters), 
               label.r = unit(0.1, 'cm'), fill='white', alpha=0.8, 
               size=5/ggplot2::.pt, label.padding=unit(0.18, 'lines')) +
    scale_x_continuous(expand=c(0.5,0.5)) +
    scale_y_continuous(expand=c(0.5,0.5)) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nect_cluster_labels_umap.pdf', width=2.85*2, height=1.9*2, unit='cm')


plot_df <- nect_meta %>% 
  mutate(class_show=case_when(
    seurat_clusters==5&UMAP_2<0 ~ 'NNE',
    seurat_clusters%in%c(1,2) ~ 'NE',
    seurat_clusters%in%c(4) ~ 'PSC',
    T ~ 'other'
  ))


nect_colors <- c('#90a4ae', '#e0e0e0', '#00897b', '#455a64')
names(nect_colors) <- c('other', 'PSC', 'NE', 'NNE')

ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=class_show)) +
  geom_point(size=3, shape=21, stroke=0.05) +
  scale_fill_manual(values=nect_colors) +
  scale_alpha_manual(values=c(0.5,1)) +
  theme_void() + no_legend()

ggsave('plots/paper/supp_pat/nect_annot_umap.png', width=15, height=10, unit='cm')




## DE heatmap 
nect_cluster_mean <- t(aggregate_matrix(t(nect_cluster_expr), nect$seurat_clusters))

vulcano_plot(nect_cluster_de, top_only = F)

nect_cluster_top <- nect_cluster_de %>% 
    filter(feature%in%good_genes, padj<1e-5) %>% 
    group_by(group) %>% arrange(group) %>% 
    top_n(4, fc) 
    
nect_cluster_marker_expr <- nect_cluster_mean[unique(nect_cluster_top$feature), ]

cluster_order <- nect_cluster_marker_expr %>% t() %>% dist() %>% hclust() %>% {.$label[.$order]}

nect_marker_expr <- nect_cluster_marker_expr %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='cluster', values_to='expr') %>% 
    inner_join(nect_cluster_top) %>% 
    mutate(group=factor(group, levels=cluster_order)) %>% 
    arrange(group) %>% 
    group_by(feature) %>% 
    mutate(
        expr=scale01(expr), 
        cluster=factor(cluster, levels=cluster_order),
        feature=factor(feature, levels=rev(unique(.$feature)))
    )

ggplot(nect_marker_expr, aes(cluster, feature, fill=expr)) +
    geom_tile() +
    scale_x_discrete(expand=c(0,0), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::brewer.greys(100)) +
    article_text() + no_legend() + 
    theme(
        axis.ticks.x = element_blank(),
        panel.border = element_blank()
    ) +
    labs(x='Louvain cluster', y='Gene')
    
ggsave('plots/paper/supp_pat/nect_marker_heatmap.pdf', width=6, height=5, unit='cm')



#### Nect NNE cluster ####
feature_plot(nect, features=c('PAX6', 'DLX6', 'TFAP2C', 'GATA3', 'HES5', 'LHX2', 'ZIC1', 'SOX2', 'SIX3', 'DLX3', 'DLX5', 'VIM',
                              'SOX3', 'POU3F1', 'GLI3', 'SIX1', 'SIX4', 'EYA1', 'EYA2', 'PAX2', 'PAX8', 'CLU'), order=T) &
  scale_color_gradientn(colors=gyylgnbu())

p1 <- feature_plot(nect, features=c('DLX5'), order=T, pt.size=0.2)
p2 <- feature_plot(nect, features=c('SOX3'), order=T, pt.size=0.2)
p3 <- feature_plot(nect, features=c('TFAP2A'), order=T, pt.size=0.2)
p4 <- feature_plot(nect, features=c('LHX2'), order=T, pt.size=0.2)

(p1 | p2 | p3 | p4) &
  scale_color_gradientn(colors=gyylgnbu()) &
  theme(plot.title=element_blank(), plot.margin = margin(0,0,0,0))

ggsave('plots/paper/supp_pat/nect_nne_marker_umap.png', width=16, height=3, unit='cm')


### Vulcano plot 
plot_df <- nect_cluster_de %>% 
  filter(group==5) %>% 
  mutate(label=feature%in%c('VIM', 'CLU', 'DLX5', 'TFAP2A', 'KRT18'))

ggplot(filter(plot_df, !label), aes(fc, -log10(padj), label=feature, alpha=padj<1e-5, color=padj<1e-5)) +
  geom_point(size=0.05, shape=16) +
  geom_point(data=filter(plot_df, label), color=red, size=0.3) +
  geom_text_repel(data=filter(plot_df, label), size=5/ggplot2::.pt) +
  scale_color_manual(values=c('grey', 'black')) +
  theme_rangeframe() + scale_axis_rangeframe() + no_legend() + article_text() +
  labs(x='Log fold change', y=expression(-log[10](FDR)))

ggsave('plots/paper/supp_pat/nect_cl5_marker_vulcano.pdf', width=3, height=3.2, unit='cm')




#### Embed early-mid patterning stage (neuepithelium) ####

rnatac$nepi <- rnatac$age > 11 & rnatac$age < 30 & rnatac$pseudotime_ranks<0.5 & rnatac$pseudotime_ranks>0.15 & rnatac$lineage!='other'

dim_plot(rnatac, group.by='nepi') +
    scale_color_manual(values=c('darkgray', 'black'))


dim_plot(rnatac, group.by='lineage')
feature_plot(rnatac, features='pseudotime_ranks')

nepi <- subset(rnatac, nepi==TRUE)

nepi@active.assay <- 'RNA'
nepi <- FindVariableFeatures(nepi, nfeatures=2000)
var_feats <- nepi@assays$RNA@meta.features %>% 
    as_tibble(rownames='gene') %>%
    filter(gene%in%good_genes) %>% 
    top_n(100, vst.variance.standardized) %>% 
    pull(gene) %>% unique()
VariableFeatures(nepi) <- var_feats %>% union(organizer_markers) %>% setdiff(cc_genes_all)
# VariableFeatures(nepi) <- organizer_markers %>% setdiff(cc_genes_all)
nepi <- CellCycleScoring(nepi, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
nepi <- ScaleData(nepi, vars.to.regress = c('S.Score', 'G2M.Score'))
nepi <- nepi %>% RunPCA() 

nepi <- cluster_sim_spectrum(
    nepi, 
    label_tag='age', 
    use_dr='pca', dims_use = 1:10,
    reduction.name = 'css', reduction.key = 'CSS_'
)

css_ <- CreateAssayObject(t(nepi[['css']]@cell.embeddings))
nepi[['css']] <- NULL
nepi[['css']] <- css_
nepi <- ScaleData(nepi, vars.to.regress = c('S.Score', 'G2M.Score'), assay='css')
css_cc <- t(nepi[['css']]@scale.data)
nepi[['css']] <- NULL
nepi[['css']] <- CreateDimReducObject(as.matrix(css_cc), key = 'CSS_')

nepi <- RunUMAP(
    object = nepi,
    spread = 0.5,
    min.dist = 0.2,
    reduction = 'css', dims = 1:ncol(Reductions(nepi, 'css')),
    reduction.name = 'umap',
    reduction.key = 'UMAP_'
)

nepi <- nepi %>% FindNeighbors(reduction = 'css', dims = 1:ncol(Reductions(nepi, 'css'))) 
nepi <- nepi %>% FindClusters(resolution=1)
nepi_cluster_expr <- nepi@assays$RNA@data
nepi_cluster_mean <- t(aggregate_matrix(t(nepi_cluster_expr), nepi$seurat_clusters))
nepi_meta <- nepi@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(Reductions(nepi, 'umap')@cell.embeddings, rownames='cell'))

dim_plot(nepi, group.by=c('age', 'line', 'lineage'))

#### Write result ####
nepi %>% write_rds('data/RNA_ATAC/subsets/RNA_nepi_srt.rds')
nepi <- read_rds('data/RNA_ATAC/subsets/RNA_nepi_srt.rds')
nepi_rna <- DietSeurat(nepi, assays='RNA', dimreducs = c('umap', 'fa2', 'css'))
nepi_rna %>% SaveH5Seurat(filename = 'data/RNA_ATAC/subsets/RNA_nepi.h5Seurat', overwrite=T)
nepi_rna <- LoadH5Seurat(file = 'data/RNA_ATAC/subsets/RNA_nepi.h5Seurat')
Convert(source = 'data/RNA_ATAC/subsets/RNA_nepi.h5Seurat', dest='h5ad')



#### Correlate with organizers ####
dim_plot(nepi, group.by=c('age', 'line', 'lineage'))

genes_use <- intersect(organizer_markers, rownames(nepi_cluster_mean))
organizer_cor <- cor(as.matrix(nepi_cluster_mean[genes_use, ]), as.matrix(organizer_mean[genes_use, ]))
pheatmap::pheatmap(organizer_cor)

organizer_cor_df <- organizer_cor %>% 
    as_tibble(rownames='seurat_clusters') %>% pivot_longer(!seurat_clusters, names_to='organizer', values_to='corr') %>% 
    inner_join(nepi_meta)


p2 <- dim_plot(nepi, group.by=c('age')) +
    scale_color_manual(values=age_colors)
p3 <- dim_plot(nepi, group.by=c('lineage')) +
    scale_color_manual(values=lineage_colors)

(p2 + p3) / p1 + plot_layout(heights=c(1,2))



#### DE between clusters ####

nepi_cluster_de <- de(nepi, groups = 'seurat_clusters')
nepi_cluster_de <- nepi_cluster_de %>% 
    filter(gene%in%good_genes)

p1 <- dim_plot(nepi, label=T)
p2 <- dim_plot(nepi, group.by=c('lineage')) +
    scale_color_manual(values=lineage_colors)
p3 <- dim_plot(nepi, group.by=c('age')) +
    scale_color_manual(values=age_colors)
p4 <- ggplot(organizer_cor_df, aes(UMAP_1, UMAP_2, color=corr)) +
    geom_point(size=0.3) +
    scale_color_gradientn(colors=gyrdpu()) +
    facet_wrap(~organizer) + theme_void()

p5 <- vulcano_plot(nepi_cluster_de)

(p1 + p2 + p3) / p4 / p5 + plot_layout(heights=c(1,3,4))
ggsave('plots/RNA_ATAC/subsets/nepi_overview.png', width=14, height=18)




feature_plot(nepi, features=c('IGFBP5', 'WNT9A', 'WNT3A', 'PTCH1', 'SIX6', 'LMX1A', 'DMRT3', 'BNC2', 'TOX2', 'MSX2', 'FEZF2', 'RAX', 'SIX3', 'PTX3'), order=T)


#### Plots for figure ####

rnatac_meta <- rnatac@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(Reductions(rnatac, 'umap')@cell.embeddings, rownames='cell'))

plot_df <- rnatac_meta
ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=nepi, alpha=nepi)) +
    geom_point(size=2, shape=21, stroke=0.1) +
    scale_fill_manual(values=c('darkgrey', '#212121')) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/full_nepi_umap.png', width=20, height=12, unit='cm')


plot_df <- nepi_meta
ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=line)) +
    geom_point(size=2, shape=21, stroke=0.1) +
    scale_fill_manual(values=line_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nepi_line_umap.png', width=16, height=12, unit='cm')

ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=factor(age))) +
    geom_point(size=2, shape=21, stroke=0.1) +
    scale_fill_manual(values=age_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nepi_age_umap.png', width=16, height=12, unit='cm')


ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=factor(lineage))) +
    geom_point(size=2, shape=21, stroke=0.1) +
    scale_fill_manual(values=lineage_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nepi_lineage_umap.png', width=16, height=12, unit='cm')

ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=pseudotime_ranks)) +
    geom_point(size=2, shape=21, stroke=0.1) +
    scale_fill_viridis(option='magma', limits=c(0,1), direction=-1) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nepi_pt_umap.png', width=16, height=12, unit='cm')

cluster_umap_means <- Reductions(nepi, 'umap')@cell.embeddings %>% 
  aggregate_matrix(nepi$seurat_clusters) %>% as_tibble(rownames='seurat_clusters')
plot_df <- nepi_meta 
ggplot(plot_df, aes(UMAP_1, UMAP_2, color=seurat_clusters)) +
  geom_point(size=1, alpha=0.5) +
  theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nepi_clusters_umap.png', width=2.85*3, height=2.2*3, unit='cm')


ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=seurat_clusters)) +
  geom_label(data=cluster_umap_means, mapping=aes(label=seurat_clusters), 
             label.r = unit(0.1, 'cm'), fill='white', alpha=0.8, 
             size=5/ggplot2::.pt, label.padding=unit(0.18, 'lines')) +
  scale_x_continuous(expand=c(0.5,0.5)) +
  scale_y_continuous(expand=c(0.5,0.5)) +
  theme_void() + no_legend()
ggsave('plots/paper/supp_pat/nepi_cluster_labels_umap.pdf', width=2.85*2, height=2.2*2, unit='cm')


## DE heatmap 
nepi_cluster_mean <- t(aggregate_matrix(t(nepi_cluster_expr), nepi$seurat_clusters))

nepi_cluster_top <- nepi_cluster_de %>% 
    filter(feature%in%good_genes, padj<1e-5) %>% 
    group_by(group) %>% arrange(group) %>% 
    top_n(3, fc) 

organizer_maerkers <- c('RAX', 'SIX3', 'SIX6',
                        'DLX2', 'GSX2', 'PTX3',
                        'LMX1A', 'MSX1', 'DMRT3', 'WNT3A')


genes_plot <- union(nepi_cluster_top$feature, organizer_maerkers)
nepi_cluster_marker_expr <- nepi_cluster_mean[genes_plot, ]

cluster_order <- nepi_cluster_marker_expr %>% t() %>% dist() %>% hclust(method='ward.D2') %>% {.$label[.$order]}
gene_order <- nepi_cluster_marker_expr %>% dist() %>% hclust(method='ward.D2') %>% {.$label[.$order]}
gene_order <- c('LDHA', 'BRF3', 'FTL', 'BTF3', 'IGFPB5', 'PTN', 'L1TD1', 'IGFBP5', 'TPBG', 'ID3', 'WNT3A', 'LMX1A', 'DMRT3', 'MSX1', 'WLS', 'PRTG', 'DLK1', 'IGFBP2', 'ODC1', 
                'VIM', 'CKB', 'HES5', 'GSX2', 'DLX2', 'PTX3', 'IGFPB2', 'FGF8', 'SIX3', 'NNAT', 'IGFBP3', 'SFRP2', 'RAX', 'SIX6', 'NR2F1', 'CRABP1') %>% rev()

nepi_marker_expr <- nepi_cluster_marker_expr %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='cluster', values_to='expr') %>% 
    # inner_join(nepi_cluster_top) %>% 
    mutate(cluster=factor(cluster, levels=cluster_order), feature=factor(feature, levels=gene_order)) %>% 
    group_by(feature) %>% 
    mutate(
        expr=scale01(expr)
    )

ggplot(nepi_marker_expr, aes(cluster, feature, fill=expr)) +
    geom_tile() +
    scale_x_discrete(expand=c(0,0), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::brewer.greys(100)) +
    article_text() + no_legend() + 
    theme(
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(x='Louvain cluster', y=NULL)

ggsave('plots/paper/supp_pat/nepi_marker_heatmap.pdf', width=6, height=7, unit='cm')


## Organizer cor heatmap

org_order <- c(
    'midbrain_floor_plate', 'hindbrain_floor_plate', 'midbrain_basal_plate', 'zona_limitans', 'hypothalamic_floor_plate',
    'roof_plate', 'telencephalic_roof_plate', 'cortical_hem', 'hindbrain_roof_plate',
    'isthmus', 'ac_pole','antihem'
)

org_colors <- c(
    rep(yellow, 5), rep(purple, 4), rep(blue, 2), red
)
names(org_colors) <- org_order

cluster_clust <- organizer_cor %>% dist() %>% hclust(method='ward.D2')
cluster_order <- cluster_clust %>% {.$label[.$order]}

plot_df <- organizer_cor_df %>% 
    distinct(seurat_clusters, organizer, corr) %>% 
    mutate(
        seurat_clusters=factor(seurat_clusters, levels=cluster_order),
        organizer=factor(organizer, levels=rev(org_order))
    ) %>% group_by(organizer) %>% 
    mutate(scale_corr=scale01(corr))
    

p1r <- ggplot(plot_df, aes(seurat_clusters, organizer, fill=corr)) +
    geom_tile() +
    scale_x_discrete(expand=c(0,0), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::brewer.bupu(100)) +
    article_text() + no_legend() +
    theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,0)
    )


p2 <- ggplot(plot_df, aes(1, organizer, fill=organizer)) +
    geom_tile() + no_legend() + article_text() +
    scale_fill_manual(values=org_colors) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    no_x_text() +
    theme(
        axis.title = element_blank(),
        plot.margin = margin(0,0,0,0)
    )

p2 + p1r + plot_layout(guides = 'collect', widths=c(1,15))

ggsave('plots/paper/supp_pat/nepi_organizer_corr_heatmap.pdf', width=6, height=3.8, unit='cm')


print_scale(pals::brewer.bupu(100))
ggsave('plots/paper/supp_pat/bupu.pdf', width=6, height=3.8, unit='cm')

print_scale(pals::brewer.rdpu(100))
ggsave('plots/paper/supp_pat/rdpu.pdf', width=6, height=3.8, unit='cm')


## Organizer cor feature plots 
organizer_cor_meta <- organizer_cor_df %>% select(cell, organizer, corr) %>% pivot_wider(names_from=organizer, values_from=corr) %>% 
    column_to_rownames('cell')
nepi <- AddMetaData(nepi, organizer_cor_meta)


p1 <- feature_plot(nepi, features='hypothalamic_floor_plate', pt.size=0.2, order=T) &
    scale_color_gradientn(colors=pals::brewer.rdpu(100)) &
    theme(
        plot.title = element_blank(),
        plot.margin = margin(0,0,0,0)
    )
p2 <- feature_plot(nepi, features='roof_plate', pt.size=0.2, order=T) &
    scale_color_gradientn(colors=pals::brewer.rdpu(100)) &
    theme(
        plot.title = element_blank(),
        plot.margin = margin(0,0,0,0)
    )
p3 <- feature_plot(nepi, features='antihem', pt.size=0.2, order=T) &
    scale_color_gradientn(colors=pals::brewer.rdpu(100)) &
    theme(
        plot.title = element_blank(),
        plot.margin = margin(0,0,0,0)
    )

p1 | p2 | p3
ggsave('plots/paper/supp_pat/nepi_organizer_corr_umap.png', width=24, height=6, unit='cm')





ggplot(organizer_cor_df, aes(UMAP_1, UMAP_2, color=corr)) +
    geom_point(size=0.3) +
    scale_color_gradientn(colors=gyrdpu()) +
    facet_wrap(~organizer) + theme_void()

ggplot(organizer_cor_df, aes(UMAP_1, UMAP_2, color=corr)) +
    geom_point(size=0.3) +
    scale_color_gradientn(colors=gyrdpu()) +
    facet_wrap(~organizer) + theme_void()





#### Embed mid branching stage (branch) ####

rnatac$branch <- rnatac$age > 11 & rnatac$pseudotime_ranks<0.85 & rnatac$pseudotime_ranks>0.15 & rnatac$lineage!='other'

dim_plot(rnatac, group.by='branch') +
    scale_color_manual(values=c('darkgray', 'black'))


dim_plot(rnatac, group.by='lineage')
feature_plot(rnatac, features='pseudotime_ranks')

branch <- subset(rnatac, branch==TRUE)

branch@active.assay <- 'RNA'
branch <- FindVariableFeatures(branch, nfeatures=2000)
var_feats <- branch@assays$RNA@meta.features %>% 
    as_tibble(rownames='gene') %>%
    filter(gene%in%good_genes) %>% 
    top_n(400, vst.variance.standardized) %>% 
    pull(gene) %>% unique()
VariableFeatures(branch) <- var_feats %>% union(organizer_markers) %>% setdiff(cc_genes_all)
# VariableFeatures(branch) <- var_feats %>% setdiff(cc_genes_all)
# VariableFeatures(branch) <- organizer_markers %>% setdiff(cc_genes_all)
branch <- CellCycleScoring(branch, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
branch <- ScaleData(branch, vars.to.regress = c('S.Score', 'G2M.Score'))
branch <- branch %>% RunPCA() 

branch <- cluster_sim_spectrum(
    branch, 
    label_tag='age', 
    use_dr='pca', dims_use = 1:10,
    reduction.name = 'css', reduction.key = 'CSS_'
)

css_ <- CreateAssayObject(t(branch[['css']]@cell.embeddings))
branch[['css']] <- NULL
branch[['css']] <- css_
branch <- ScaleData(branch, vars.to.regress = c('S.Score', 'G2M.Score'), assay='css')
css_cc <- t(branch[['css']]@scale.data)
branch[['css']] <- NULL
branch[['css']] <- CreateDimReducObject(as.matrix(css_cc), key = 'CSS_')

branch <- RunUMAP(
    object = branch,
    spread = 0.5,
    min.dist = 0.2,
    reduction = 'css', dims = 1:ncol(Reductions(branch, 'css')),
    reduction.name = 'umap',
    reduction.key = 'UMAP_'
)

branch <- branch %>% FindNeighbors(reduction = 'css', dims = 1:ncol(Reductions(branch, 'css'))) 
branch <- branch %>% FindClusters(resolution=1)
branch_cluster_expr <- branch@assays$RNA@data
branch_cluster_mean <- t(aggregate_matrix(t(branch_cluster_expr), branch$seurat_clusters))
branch_meta <- branch@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(Reductions(branch, 'umap')@cell.embeddings, rownames='cell'))

dim_plot(branch, group.by=c('age', 'line', 'lineage'))



branch %>% write_rds('data/RNA_ATAC/subsets/RNA_branch.srt')
branch_rna <- DietSeurat(branch, assays='RNA', dimreducs = c('umap', 'fa2', 'css'))
branch_rna %>% SaveH5Seurat(filename = 'data/RNA_ATAC/subsets/RNA_branch.h5Seurat', overwrite=T)
Convert(source = 'data/RNA_ATAC/subsets/RNA_branch.h5Seurat', dest='h5ad')



#### Correlate with organizers ####
organizer_markers <- organizer_de %>% 
    filter(padj<0.01) %>% group_by(group) %>% top_n(30, fc) %>% pull(gene)

genes_use <- intersect(organizer_markers, rownames(branch_cluster_mean))
organizer_cor <- cor(as.matrix(branch_cluster_mean[genes_use, ]), as.matrix(organizer_mean[genes_use, ]))
pheatmap::pheatmap(organizer_cor)

organizer_cor_df <- organizer_cor %>% 
    as_tibble(rownames='seurat_clusters') %>% pivot_longer(!seurat_clusters, names_to='organizer', values_to='corr') %>% 
    inner_join(branch_meta)


p2 <- dim_plot(branch, group.by=c('age')) +
    scale_color_manual(values=age_colors)
p3 <- dim_plot(branch, group.by=c('lineage')) +
    scale_color_manual(values=lineage_colors)

(p2 + p3) / p1 + plot_layout(heights=c(1,2))



#### DE between clusters ####

branch_cluster_de <- de(branch, groups = 'seurat_clusters')
branch_cluster_de <- branch_cluster_de %>% 
    filter(gene%in%good_genes)

p1 <- dim_plot(branch, label=T)
p2 <- dim_plot(branch, group.by=c('lineage')) +
    scale_color_manual(values=lineage_colors)
p3 <- dim_plot(branch, group.by=c('age')) +
    scale_color_manual(values=age_colors)
p3.1 <- dim_plot(branch, group.by=c('line')) +
    scale_color_manual(values=line_colors)
p4 <- ggplot(organizer_cor_df, aes(UMAP_1, UMAP_2, color=corr)) +
    geom_point(size=0.3) +
    scale_color_gradientn(colors=gyrdpu()) +
    facet_wrap(~organizer) + theme_void()

p5 <- vulcano_plot(branch_cluster_de)

(p1 | p2 | p3 | p3.1) / p4 / p5 + plot_layout(heights=c(1,3,4))
ggsave('plots/RNA_ATAC/subsets/branch_overview.png', width=14, height=18)



#### Plots for figure ####
rnatac_meta <- rnatac@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(Reductions(rnatac, 'umap')@cell.embeddings, rownames='cell'))

plot_df <- rnatac_meta
ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=branch, alpha=branch)) +
    geom_point(size=2, shape=21, stroke=0.1) +
    scale_fill_manual(values=c('darkgrey', '#212121')) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/full_branch_umap.png', width=20, height=12, unit='cm')


plot_df <- sample_n(branch_meta, nrow(branch_meta))
ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=line)) +
    geom_point(size=1.5, shape=21, stroke=0.1) +
    scale_fill_manual(values=line_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/branch_line_umap.png', width=16, height=12, unit='cm')

ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=factor(age))) +
    geom_point(size=1.5, shape=21, stroke=0.1) +
    scale_fill_manual(values=age_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/branch_age_umap.png', width=16, height=12, unit='cm')

plot_df <- sample_n(branch_meta, nrow(branch_meta)) %>% arrange(lineage%in%c('telencephalon', 'nt', 'ctx', 'early'))
ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=factor(lineage))) +
    geom_point(size=2, shape=21, stroke=0.1) +
    scale_fill_manual(values=lineage_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/branch_lineage_umap.png', width=16, height=12, unit='cm')

ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=pseudotime_ranks)) +
    geom_point(size=2, shape=21, stroke=0.1, alpha=0.5) +
    scale_fill_viridis(option='magma', limits=c(0,1), direction=-1) +
    theme_void() + no_legend()
ggsave('plots/paper/supp_pat/branch_pt_umap.png', width=16, height=12, unit='cm')

plot_df <- branch_meta
p1 <- ggplot(arrange(plot_df, to_nt_ranks), aes(UMAP_1, UMAP_2, color=to_nt_ranks)) +
    geom_point(size=1, stroke=0.05) +
    scale_color_gradientn(colors=rev(pals::ocean.deep(100))) +
    theme_void() + no_legend()

p2 <- ggplot(arrange(plot_df, to_ge_ranks), aes(UMAP_1, UMAP_2, color=to_ge_ranks)) +
    geom_point(size=1, stroke=0.05) +
    scale_color_gradientn(colors=rev(pals::ocean.deep(100))) +
    theme_void() + no_legend()

p3 <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color=to_ctx_ranks)) +
    geom_point(size=1, stroke=0.05) +
    scale_color_gradientn(colors=rev(pals::ocean.deep(100))) +
    theme_void() + no_legend()

p1 + p2 + p3
ggsave('plots/paper/supp_pat/branch_trans_prob_umap.png', width=30, height=8, unit='cm')



cluster_umap_means <- Reductions(branch, 'umap')@cell.embeddings %>% 
  aggregate_matrix(branch$seurat_clusters) %>% as_tibble(rownames='seurat_clusters')
plot_df <- branch_meta 
ggplot(plot_df, aes(UMAP_1, UMAP_2, color=seurat_clusters)) +
  geom_point(size=0.8, alpha=0.5) +
  theme_void() + no_legend()
ggsave('plots/paper/supp_pat/branch_clusters_umap.png', width=2.4*3, height=1.8*3, unit='cm')


ggplot(plot_df, aes(UMAP_1, UMAP_2, fill=seurat_clusters)) +
  geom_label(data=cluster_umap_means, mapping=aes(label=seurat_clusters), 
             label.r = unit(0.1, 'cm'), fill='white', alpha=0.8, 
             size=5/ggplot2::.pt, label.padding=unit(0.18, 'lines')) +
  scale_x_continuous(expand=c(0.5,0.5)) +
  scale_y_continuous(expand=c(0.5,0.5)) +
  theme_void() + no_legend()
ggsave('plots/paper/supp_pat/branch_cluster_labels_umap.pdf', width=2.4*2, height=1.8*2, unit='cm')


## DE heatmap 
branch_cluster_mean <- t(aggregate_matrix(t(branch_cluster_expr), branch$seurat_clusters))

branch_cluster_top <- branch_cluster_de %>% 
    group_by(feature) %>% filter(fc==max(fc)) %>% 
    filter(feature%in%good_genes, padj<1e-5) %>% 
    group_by(group) %>% 
    arrange(group, desc(fc)) %>% 
    filter(row_number()<=2 | feature%in%c('BMP7', 'FGF8', 'WNT8B'))

branch_cluster_marker_expr <- branch_cluster_mean[unique(branch_cluster_top$feature), ]

cluster_order <- branch_cluster_marker_expr %>% t() %>% dist() %>% hclust() %>% {.$label[.$order]}

branch_marker_expr <- branch_cluster_marker_expr %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='cluster', values_to='expr') %>% 
    inner_join(branch_cluster_top) %>% 
    mutate(group=factor(group, levels=cluster_order)) %>% 
    arrange(group) %>% 
    group_by(feature) %>% 
    mutate(
        expr=scale01(expr), 
        cluster=factor(cluster, levels=cluster_order),
        feature=factor(feature, levels=rev(unique(.$feature)))
    )

ggplot(branch_marker_expr, aes(cluster, feature, fill=expr)) +
    geom_tile() +
    scale_x_discrete(expand=c(0,0), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::brewer.greys(100)) +
    article_text() + no_legend() + 
    theme(
        axis.ticks.x = element_blank(),
        panel.border = element_blank()
    ) +
    labs(x='Louvain cluster', y='Gene')

ggsave('plots/paper/supp_pat/branch_marker_heatmap.pdf', width=6, height=5.5, unit='cm')



