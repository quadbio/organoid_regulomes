require(plyr, quietly=T)
require(tidyverse, quietly=T)
require(Seurat, quietly=T)
require(Matrix, quietly=T)
require(umap, quietly=T)
require(uwot, quietly=T)
require(Rtsne, quietly=T)
require(phateR, quietly=T)
require(irlba, quietly=T)
require(gridExtra, quietly=T)


basic_qc <- function(seurat_object, file=NULL){

    # Basic QC
    mito_genes <- grep(
        pattern = "^MT-",
        x = rownames(seurat_object@assays$RNA@counts), value=T)
    mito_expr <- seurat_object@assays$RNA@counts[mito_genes, ]
    percent_mito <- colSums(mito_expr)/colSums(seurat_object@assays$RNA@counts)
    seurat_object <- AddMetaData(
        object = seurat_object,
        metadata = percent_mito,
        col.name = "percent_mito")

    p1 <- ggplot(seurat_object@meta.data, aes(nCount_RNA)) +
        geom_histogram(color="black")
    p2 <- ggplot(seurat_object@meta.data, aes(nFeature_RNA)) +
        geom_histogram(color="black")
    p3 <- ggplot(seurat_object@meta.data, aes(percent_mito)) +
        geom_histogram(color="black")

    if (!is.null(file)){
        pdf(file)
        grid.arrange(p1,p2,p3)
        dev.off()
    } else {
        grid.arrange(p1,p2,p3)
    }

    return(seurat_object)
}

umap <- function(
    expr_mat,
    n_pcs = 15,
    do_pca = T,
    method = 'uwot',
    config = umap.defaults,
    ...
){
    if (do_pca){
        pca_mat <- prcomp_irlba(expr_mat, n=n_pcs)$x
        rownames(pca_mat) <- rownames(expr_mat)
    } else {
        umap_tbl <- umap::umap(pca_mat, config=config, method=method, ...)$layout %>%
            as_tibble(rownames="cell")
    }
    if (method == 'uwot'){
        umap_tbl <- uwot::umap(as.matrix(pca_mat), ...) %>%
            {rownames(.) <- rownames(pca_mat); .} %>%
            as_tibble(rownames='cell')
    } else {
        umap_tbl <- umap::umap(pca_mat, config=config, method=method, ...)$layout %>%
            as_tibble(rownames="cell")
    }
    colnames(umap_tbl) <- c('cell', paste0('UMAP', seq(ncol(umap_tbl)-1)))
    return(umap_tbl)
}

tsne <- function(expr_mat, n_pcs=15, do_pca=T, perplexity=30){
    if (do_pca){
        pca_mat <- prcomp_irlba(expr_mat, n=n_pcs)$x
        rownames(pca_mat) <- rownames(expr_mat)
    } else {
        pca_mat <- expr_mat
    }
    tsne_tbl <- Rtsne(pca_mat, check_duplicates=F, perplexity=perplexity)$Y %>%
        as_tibble() %>%
        mutate(cell=rownames(pca_mat))
    colnames(tsne_tbl) <- c(paste0('tSNE', seq(ncol(tsne_tbl)-1)), 'cell')
    return(tsne_tbl)
}


cluster <- function(expr_mat,
                    n_pcs = 15,
                    do_pca = T,
                    method = 'leiden',
                    resolution = 0.8){
    if (do_pca){
        pca_mat <- prcomp_irlba(expr_mat, n=n_pcs)$x
        rownames(pca_mat) <- rownames(expr_mat)
    } else {
        pca_mat <- expr_mat
    }
    if (method == 'leiden'){
        alg_num <- 4
    } else if (method == 'louvain'){
        alg_num <- 1
    }
    knn <- FindNeighbors(pca_mat, verbose=F)
    clusters <- FindClusters(knn$snn, verbose=F, method=alg_num, resolution=resolution) %>%
        as_tibble(rownames='cell')
    colnames(clusters) <- c('cell', paste0(method, '_', as.character(resolution)))

    return(clusters)
}


variable_features <- function(object, ...){
    UseMethod(generic = 'variable_features', object = object)
}

variable_features.default <- function(object, method = 'vst', nfeatures = 2000){

    var_feat <- FindVariableFeatures(t(object), verbose=F) %>%
        as_tibble(rownames='gene') %>%
        arrange(desc(vst.variance.standardized)) %>%
        dplyr::mutate(highvar = dplyr::row_number() <= nfeatures) %>%
        {colnames(.) <- str_replace_all(colnames(.), '\\.', '_'); .}

    return(var_feat)
}

variable_features.Seurat <- function(
    object, method='vst', nfeatures=2000, assay='RNA', slot='counts'
){
    expr_mat <- GetAssayData(object, assay=assay, slot=slot)
    var_feat <- FindVariableFeatures(expr_mat, verbose=F) %>%
        as_tibble(rownames='gene') %>%
        arrange(desc(vst.variance.standardized)) %>%
        dplyr::mutate(highvar = dplyr::row_number() <= nfeatures) %>%
        {colnames(.) <- str_replace_all(colnames(.), '\\.', '_'); .}
    return(var_feat)
}

log_normalize <- function(expr_mat, scale_factor = 10000){
    out_mat <- t(LogNormalize(t(expr_mat), scale.factor = scale_factor, verbose = F))
    return(out_mat)
}

elbow_plot <- function(seurat_object){
    pc_df <- tibble(
        Stdev=seurat_object@reductions$pca@stdev,
        PC=seq_along(seurat_object@reductions$pca@stdev)
    )
    p <- ggplot(pc_df, aes(PC, Stdev)) +
        geom_point() +
        geom_vline(xintercept=20, color='red')
    return(p)
}

pca <- function(expr_mat, n=20, to_df=F){
    pca_x <- irlba::prcomp_irlba(expr_mat, n=n)$x
    rownames(pca_x) <- rownames(expr_mat)
    if (to_df){
        pca_x <- as_tibble(pca_x, rownames='cell')
    }
    return(pca_x)
}
