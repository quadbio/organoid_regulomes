source('~/scripts/single_cell/de.R')

setwd('~/projects/early/')

library(igraph)
library(simspec)
library(ggraph)
library(tidygraph)
library(destiny)
library(DDRTree)


#### Functions ####
knn <- function(dst, k){
    k <- k+1
    knn <- RANN::nn2(as.matrix(dst), k=k)
    knn_tbl <- knn$nn.idx[, 2:k] %>% 
        as_tibble(rownames='source') %>% 
        pivot_longer(!source, names_to='id', values_to='idx') %>% 
        inner_join(pivot_longer(as_tibble(knn$nn.dists, rownames='source'),
                                !source, names_to='id', values_to='value')) %>% 
        dplyr::select(source, idx, value) %>%
        dplyr::mutate(target=rownames(as.matrix(dst))[as.numeric(.$idx)]) %>%
        dplyr::mutate(source=rownames(as.matrix(dst))[as.numeric(.$source)]) %>%
        dplyr::select(source, target, value)
    return(knn_tbl)
}


#### Read stuff ####
rnatac <- read_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2.1_srt.rds')
# rnatac %>% write_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2.1_srt.rds')

rnatac <- FindNeighbors(rnatac, reduction='css', dims = 1:ncol(Reductions(rnatac, 'css')))
rnatac <- FindClusters(rnatac, resolution = 20)
dim_plot(rnatac, label=T, cols=many_more) + no_legend()

mid_pat <- read_rds('data/RNA_ATAC/subsets/RNA_mid_pat.srt')

cr_meta <- read_tsv('data/RNA_ATAC/velocity/RNA_ATAC_full_cellrank_probs.tsv')
colnames(cr_meta)[1] <- 'cell'
cr_meta_df <- column_to_rownames(cr_meta, 'cell')

rnatac <- AddMetaData(rnatac, cr_meta_df)
mid_pat <- AddMetaData(mid_pat, cr_meta_df)

cr_trans_mat <- readMM('data/RNA_ATAC/velocity/velo_cr_transition.mtx')
colnames(cr_trans_mat) <- colnames(rnatac)
rownames(cr_trans_mat) <- colnames(rnatac)

paga_conn_mat <- readMM('data/RNA_ATAC/velocity/velo_paga_connectivities.mtx')
colnames(paga_conn_mat) <- seq(ncol(paga_conn_mat)) - 1
rownames(paga_conn_mat) <- seq(nrow(paga_conn_mat)) - 1



#### Try dimreducs on cr space ####
cr_space <- cr_meta %>% 
    select(cell, to_ge_ranks, to_ctx_ranks, to_nt_ranks, pseudotime_ranks) %>% 
    column_to_rownames('cell') %>% as.matrix() 

cr_pca <- pca(cr_space[,1:3], n = 2, to_df = T) %>% 
    inner_join(cr_meta)

ggplot(cr_pca, aes(PC1, pseudotime_ranks, color=to_ctx_ranks)) +
    geom_point(size=0.5, alpha=0.8) +
    scale_color_gradientn(colors=gyylgnbu()) +
    theme_void()


#### Reduce to clusters ####
feature_plot(rnatac, features=c('to_ctx_ranks', 'to_ge_ranks', 'to_nt_ranks'))
dim_plot(rnatac, group.by=c('RNA_snn_res.20', 'RNA_snn_res.10', 'RNA_snn_res.2')) &
    scale_color_manual(values=many_more) & no_legend()

cluster_meta <- rnatac@meta.data[, c('RNA_snn_res.20', 'RNA_snn_res.10', 'RNA_snn_res.2')]

cluster_counts <- as.numeric(table(rnatac$RNA_snn_res.20))
cr_hc_space <- aggregate.Matrix(column_to_rownames(cr_pca, 'cell'), groupings = rnatac$RNA_snn_res.20) / cluster_counts
cr_hc_df <- cr_hc_space %>% as_tibble(rownames='cluster') %>% 
    filter(!cluster%in%c(142,162))


cluster_counts <- table(rnatac$RNA_snn_res.20, rnatac$age)[cr_hc_df$cluster, ]

age_clusters <- map(as.character(sort(unique(rnatac$age))), function(x){
    rownames(cluster_counts)[cluster_counts[, x]>0]  
})


ggplot(cr_hc_df, aes(x=.panel_x, y=.panel_y, color=to_ctx_ranks)) +
    geom_point(alpha=0.8, size=3) +
    scale_color_viridis() +
    facet_matrix(vars(to_ge_ranks, to_ctx_ranks, to_nt_ranks, pseudotime_ranks), layer.diag=FALSE) +
    theme_article()

# get_umap_coords for clusters
cr_hc_space <- aggregate.Matrix(Reductions(rnatac, 'umap')@cell.embeddings, groupings = rnatac$RNA_snn_res.20) / cluster_counts
cr_hc_space_df <- as_tibble(cr_hc_space, rownames='cluster') %>% mutate(cluster=as.numeric(cluster))
ggplot(cr_hc_space_df, aes(UMAP_1, UMAP_2)) +
    geom_point()

# Get transition probs for clusters 
cr_trans_mat <- cr_trans_mat[colnames(rnatac), colnames(rnatac)]
hc_trans_mat <- aggregate.Matrix(cr_trans_mat, groupings = rnatac$RNA_snn_res.20) / cluster_counts 
hc_trans_mat <- t(aggregate.Matrix(t(hc_trans_mat), groupings = rnatac$RNA_snn_res.20) / cluster_counts)
diag(hc_trans_mat) <- 0

cr_hc_df <- inner_join(cr_hc_df, cr_hc_space_df)
# cr_hc_df %>% write_tsv('data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv')
# cr_hc_df <- read_tsv('data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv')



#### Get NN graph ####
k <- 31
cr_hc_df$cluster <- as.character(cr_hc_df$cluster)
pt_order_idx <- order(cr_hc_df$pseudotime_ranks)
pt_order <- as.character(cr_hc_df$cluster[pt_order_idx])
cr_hc_mat <- column_to_rownames(cr_hc_df, 'cluster')[, c('to_ctx_ranks', 'to_ge_ranks', 'to_nt_ranks')]
cr_hc_mat <- cr_hc_mat[pt_order, ]

cr_hc_nn <- knn(dist(cr_hc_mat), k=k)
cr_dist <- as_tbl_graph(cr_hc_nn) %E>% as_tibble()

#### Prune graph ####
nn_mat <- as_adjacency_matrix(as_tbl_graph(cr_hc_nn))[pt_order, pt_order]

# Edges back in pseudotime
nn_mat <- triu(nn_mat)
pheatmap::pheatmap(nn_mat, cluster_cols = F, cluster_rows = F)

# By PAGA connectivities
paga_thresh <- paga_conn_mat > 0.2
pt_order %in% colnames(paga_thresh)
nn_mat <- nn_mat * paga_thresh[pt_order, pt_order]
pheatmap::pheatmap(nn_mat, cluster_cols = F, cluster_rows = F)

nn_graph <- as_tbl_graph(graph_from_adjacency_matrix(nn_mat)) %>% 
    inner_join(cr_hc_df, by=c('name'='cluster')) %>% 
    convert(to_simple) %E>%
    inner_join(cr_dist) %N>% 
    mutate(
        central_pr=centrality_pagerank() 
    )

#### Plot graph ####
ggraph(nn_graph, x=pseudotime_ranks, y=PC1) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=pseudotime_ranks, size=central_pr)) +
    geom_node_text(aes(label=name)) +
    scale_color_viridis(option='magma', direction=-1)

ggraph(nn_graph, layout='fr') +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=to_ge_ranks), size=2) +
    scale_color_viridis()


ggplot(cr_hc_df, aes(to_ctx_ranks, pseudotime_ranks, color=cluster==56)) +
    geom_point(size=0.5, alpha=0.8) +
    theme_void()

pheatmap::pheatmap(hc_trans_mat[pt_order, pt_order], cluster_cols = F, cluster_rows = F)


#### Define root and tips ####
nn_graph <- read_rds('data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')
cr_hc_df <- read_tsv('data/RNA_ATAC/lineages/res_20_clusters_cellrank_meta.tsv')
pt_order_idx <- order(cr_hc_df$pseudotime_ranks)
pt_order <- as.character(cr_hc_df$cluster[pt_order_idx])

root_cell <- cr_hc_df$cluster[which.min(cr_hc_df$pseudotime_ranks)]
ctx_tip <- filter(cr_hc_df, to_ctx_ranks>0.95) %>% {.$cluster[which.max(.$pseudotime_ranks)]} %>% as.character()
ge_tip <- filter(cr_hc_df, to_ge_ranks>0.95) %>% {.$cluster[which.max(.$pseudotime_ranks)]} %>% as.character()
nt_tip <- filter(cr_hc_df, to_nt_ranks>0.95) %>% {.$cluster[which.max(.$pseudotime_ranks)]} %>% as.character()

rwalks <- map(list(nt=nt_tip, ctx=ctx_tip, ge=ge_tip), function(tip){
    map(1:10000, function(.){
        names(random_walk(nn_graph, start=tip, steps=200, stuck='return', mode='in'))
    })
})

freqs <- map(rwalks, function(w){
    fq <- w %>% 
        map(function(x) t(as.matrix(table(x))))
    fq <- Matrix(do.call(rbind.fill.matrix, fq), sparse=T)
    fq[is.na(fq)] <- 0
    return(t(colSums(fq)) / 10000)
})
rwalk_freqs <- Matrix(do.call(rbind.fill.matrix, freqs), sparse=T)
rwalk_freqs[is.na(rwalk_freqs)] <- 0
rownames(rwalk_freqs) <- names(freqs)

reached_points <- pt_order[pt_order%in%colnames(rwalk_freqs)]
rwalk_freqs <- rwalk_freqs[, reached_points]

pheatmap::pheatmap(rwalk_freqs)

rnatac$this <- rnatac$RNA_snn_res.20==root_cell
dim_plot(rnatac, label=T, group.by='this') + no_legend()

m <- 0.01
assignments <- apply(rwalk_freqs, 2, function(i){
    if (length(i[i!=0])==1){
        return(names(i)[i!=0])
    }
    combs <- combn(i, 2, simplify = F)
    ratios <- map_dbl(combs, function(x) min(x)/max(x))
    if (all(ratios>(m/2))){
        return('early')
    }
    if (all(ratios<m)){
        return(names(which.max(i)))
    }
    high_idx <- which(ratios>m)
    high_combs <- combs[high_idx]
    high_lins <- purrr::reduce(map(high_combs, function(x) names(x)), unique)
    return(paste0(high_lins, collapse='_'))
})
assignments[assignments%in%c('ctx_ge', 'ge_ctx')] <- 'telencephalon'


# Assign unreched points 
unreached_points <- pt_order[!pt_order%in%colnames(rwalk_freqs)]

unreached_assign <- map_chr(unreached_points, function(x){
    nn <- names(neighbors(nn_graph, x, mode = 'in'))
    assign <- assignments[names(sort(paga_conn_mat[x, nn], decreasing=T)[1])]
    return(assign)
})

names(unreached_assign) <- unreached_points

assign_meta <- enframe(c(assignments, unreached_assign), 'name', 'lineage')

nn_assign_graph <- nn_graph %N>% 
    select(-lineage) %>%
    left_join(assign_meta)

set.seed(111)
ggraph(nn_assign_graph, layout='drl') +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=lineage), size=2) +
    scale_color_manual(values=lineage_colors)

set.seed(111)
ggraph(nn_assign_graph, layout='fr') +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=lineage), size=2) +
    scale_color_manual(values=lineage_colors)

ggraph(nn_assign_graph, x=pseudotime_ranks, y=PC1) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=lineage), size=3) +
    geom_node_text(aes(label=name), size=4) +
    scale_color_manual(values=lineage_colors) +
    theme_void()

nn_assign_graph %>% write_rds('data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')


#### Turn into tree ####
# Draw individual lineages
early_branch_nodes <- nn_assign_graph %>% 
    filter(lineage=='early') %N>% 
    as_tibble() %>% arrange(pseudotime_ranks) %>% mutate(pseudotime_align=pseudotime_ranks)
early_branch <- create_path(nrow(early_branch_nodes), directed = TRUE) %>% 
    mutate(name=early_branch_nodes$name) %>% inner_join(early_branch_nodes)

nt_branch_nodes <- nn_assign_graph %>% 
    filter(lineage=='nt') %N>% 
    as_tibble() %>% arrange(pseudotime_ranks) %>% mutate(pseudotime_align=pseudotime_ranks+(max(early_branch_nodes$pseudotime_ranks)-min(pseudotime_ranks)))
nt_branch <- create_path(nrow(nt_branch_nodes), directed = TRUE) %>% 
    mutate(name=nt_branch_nodes$name) %>% inner_join(nt_branch_nodes)

telencephalon_branch_nodes <- nn_assign_graph %>% 
    filter(lineage=='telencephalon') %N>% 
    as_tibble() %>% arrange(pseudotime_ranks) %>% mutate(pseudotime_align=pseudotime_ranks+(max(early_branch_nodes$pseudotime_ranks)-min(pseudotime_ranks)))
telencephalon_branch <- create_path(nrow(telencephalon_branch_nodes), directed = TRUE) %>% 
    mutate(name=telencephalon_branch_nodes$name) %>% inner_join(telencephalon_branch_nodes)

ge_branch_nodes <- nn_assign_graph %>% 
    filter(lineage=='ge') %N>% 
    as_tibble() %>% arrange(pseudotime_ranks) %>% mutate(pseudotime_align=pseudotime_ranks+(max(telencephalon_branch_nodes$pseudotime_align)-min(pseudotime_ranks)))
ge_branch <- create_path(nrow(ge_branch_nodes), directed = TRUE) %>% 
    mutate(name=ge_branch_nodes$name) %>% inner_join(ge_branch_nodes)

ctx_branch_nodes <- nn_assign_graph %>% 
    filter(lineage=='ctx') %N>% 
    as_tibble() %>% arrange(pseudotime_ranks) %>% mutate(pseudotime_align=pseudotime_ranks+(max(telencephalon_branch_nodes$pseudotime_align)-min(pseudotime_ranks)))
ctx_branch <- create_path(nrow(ctx_branch_nodes), directed = TRUE) %>% 
    mutate(name=ctx_branch_nodes$name) %>% inner_join(ctx_branch_nodes)

# Join lineages
combined_tree <- bind_graphs(early_branch, telencephalon_branch, nt_branch, ge_branch, ctx_branch) %>% 
    bind_edges(data.frame(from = early_branch_nodes$name[nrow(early_branch_nodes)], to = telencephalon_branch_nodes$name[1])) %>% 
    bind_edges(data.frame(from = early_branch_nodes$name[nrow(early_branch_nodes)], to = nt_branch_nodes$name[1])) %>% 
    bind_edges(data.frame(from = telencephalon_branch_nodes$name[nrow(telencephalon_branch_nodes)], to = ge_branch_nodes$name[1])) %>% 
    bind_edges(data.frame(from = telencephalon_branch_nodes$name[nrow(telencephalon_branch_nodes)], to = ctx_branch_nodes$name[1]))

layout_coords_x <- layout_as_tree(combined_tree)[,1]

combined_tree <- combined_tree %>% 
    mutate(coords_x = layout_coords_x)

ggraph(combined_tree, x=pseudotime_align, y=coords_x) +
    geom_edge_diagonal(alpha=0.6) +
    geom_node_point(aes(color=lineage), size=3) +
    geom_node_text(aes(label=name), size=3) +
    scale_color_manual(values=lineage_colors)

ggraph(combined_tree, x=pseudotime_align, y=coords_x) +
    geom_edge_diagonal(alpha=0.6) +
    geom_node_point(aes(color=pseudotime_ranks), size=3) +
    scale_color_viridis(option='magma', direction=-1)

ggraph(combined_tree, x=pseudotime_ranks, y=coords_x) +
    geom_edge_diagonal(alpha=0.6) +
    geom_node_point(aes(fill=to_nt_ranks), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=pals::ocean.deep(100)) +
    theme_void()
    
ggraph(combined_tree, x=pseudotime_align, y=coords_x) +
    geom_edge_diagonal(alpha=0.6) +
    geom_node_point(aes(fill=to_ctx_ranks), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=pals::ocean.deep(100)) +
    theme_void()
    
ggraph(combined_tree, x=pseudotime_align, y=coords_x) +
    geom_edge_diagonal(alpha=0.6) +
    geom_node_point(aes(fill=to_ge_ranks), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=pals::ocean.deep(100)) +
    theme_void()
    
combined_tree %>% write_rds('data/RNA_ATAC/lineages/cellrank_lineage_tree.rds')

tree_meta <- combined_tree %N>% 
    as_tibble() %>% select(name, pseudotime_align, coords_x)

nn_assign_graph <- nn_assign_graph %N>% 
    left_join(tree_meta)

nn_assign_graph %>% write_rds('data/RNA_ATAC/lineages/cellrank_lineage_graph.rds')



##### Tree layout to initiate force-directed ####
tree_coords <- nn_assign_graph %>% select(pseudotime_ranks, coords_x) %N>% as_tibble() %>% as.matrix()

set.seed(111)
p1 <- ggraph(nn_assign_graph, layout='fr') +
    geom_edge_link(alpha=0.6) +
    geom_node_point(aes(fill=lineage), size=4, shape=21, stroke=0.2) +
    scale_fill_manual(values=lineage_colors) +
    theme_void() +
    ggtitle('A - Force directed')

fr_layout <- layout_with_fr(nn_assign_graph, coords=tree_coords, niter=50) %>% 
    as_tibble() %>% mutate(name=names(V(nn_assign_graph))) %>% rename('FR1'='V1', 'FR2'='V2')

nn_assign_graph <- nn_assign_graph %N>% 
    inner_join(fr_layout)

ggraph(plot_graph, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.6) +
    geom_node_point(aes(fill=lineage), size=4, shape=21, stroke=0.2) +
    scale_fill_manual(values=lineage_colors) +
    theme_void() +
    ggtitle('B - Force directed + initiation with tree')

set.seed(111)
p2 <- ggraph(nn_assign_graph, layout='fr', coords=tree_coords, niter=50) +
    geom_edge_link(alpha=0.6) +
    geom_node_point(aes(fill=lineage), size=4, shape=21, stroke=0.2) +
    scale_fill_manual(values=lineage_colors) +
    theme_void() +
    ggtitle('B - Force directed + initiation with tree')

p3 <- ggraph(nn_assign_graph, x=pseudotime_ranks, y=PC1) +
    geom_edge_diagonal(alpha=0.6, strength=0.3) +
    geom_node_point(aes(fill=lineage), size=4, shape=21, stroke=0.2) +
    scale_fill_manual(values=lineage_colors) +
    theme_void() +
    ggtitle('C - pseudotime vs PC1')

p4 <- ggraph(nn_assign_graph, x=pseudotime_ranks, y=coords_x) +
    geom_edge_diagonal(alpha=0.6, strength=0.5) +
    geom_node_point(aes(fill=lineage), size=4, shape=21, stroke=0.2) +
    scale_fill_manual(values=lineage_colors) +
    theme_void() +
    ggtitle('D - Pseudotime tree with all edges')

p5 <- ggraph(combined_tree, x=pseudotime_align, y=coords_x) +
    geom_edge_diagonal(alpha=0.6, strength=0.5) +
    geom_node_point(aes(fill=lineage), size=4, shape=21, stroke=0.2) +
    scale_fill_manual(values=lineage_colors) +
    theme_void() +
    ggtitle('E - Aligned pseudotime tree with pruned edges')


p1 / p2 / p3 / p4 / p5 & no_legend()
ggsave('plots/RNA_ATAC/lineages/tree_graph_options.png', width=6, height=15)




















