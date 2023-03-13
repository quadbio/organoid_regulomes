require(plyr)
require(tidyverse)
require(glmnet)
require(Matrix)
require(Matrix.utils)
require(simspec)
require(irlba)
require(voxhunt)

data('brainspan')

ens_sym <- read_tsv('~/resources/mappings/ens2sym.tsv')

#### Nowakowski cell type prediction ####

if (!exists('nowa_data')){
    nowa_data <- readRDS('/links/groups/treutlein/PUBLIC_DATA/published/single_cell/2017_nowakowski_human_cortex_fetal_science/data.expr_meta_all.rds')
    nowa_meta <- as_tibble(nowa_data$meta, rownames='sample')
}

nowakowski_fit <- function(genes){
    genes_use <- intersect(rownames(nowa_data$expr), genes)
    nowa_expr <- nowa_data$expr[, nowa_meta$sample]
    nowa_expr_ranks <- t(apply(nowa_expr, 2, rank))
    fit <- cv.glmnet(nowa_expr_ranks, nowa_meta$MainCellType, family='multinomial', alpha=1)
    return(list(fit=fit, genes=genes))
}

nowakowski_predict <- function(object, newmtx){
    expr <- as.matrix(newmtx[, object$genes])
    expr_ranks <- t(apply(expr, 1, rank))
    pred <- predict(object$fit, expr_ranks, type='response')
    pred_df <- pred[,,1] %>%
        as_tibble(rownames='cell') %>%
        gather(cell_type, p_pred, -cell)
    return(pred_df)
}


#### RSS ####
rss <- function(object, ...){
    UseMethod(generic = 'rss', object = object)
}

rss.default <- function(object, only_young_samples=FALSE, do_pca=TRUE, n_pcs=20){
    if (only_young_samples){
        b10pcw <- brainspan$row_meta$fetal_before10pcw
        rss_mat <- ref_sim_spectrum(as.matrix(object), t(brainspan$matrix)[, b10pcw])
    } else {
        rss_mat <- ref_sim_spectrum(as.matrix(object), t(brainspan$matrix))
    }
    rownames(rss_mat) <- colnames(object)
    if (do_pca){
        rss_pca_mat <- prcomp_irlba(rss_mat, n=n_pcs)$x
        rownames(rss_pca_mat) <- rownames(rss_mat)
        return(rss_pca_mat)
    } else {
        return(rss_mat)
    }
}

rss.Seurat <- function(
    object,
    assay = 'RNA',
    slot = 'data',
    only_young_samples = FALSE,
    do_pca = TRUE,
    n_pcs = 20
){
    expr_mat <- GetAssayData(object, assay=assay, slot=slot)
    genes_use <- VariableFeatures(object, assay=assay)
    rss_mat <- rss(
        expr_mat[genes_use, ],
        only_young_samples = only_young_samples,
        do_pca = do_pca,
        n_pcs = n_pcs
    )
    object[['rss']] <- CreateDimReducObject(rss_mat, key='RSS_', assay=assay)
    return(object)
}
