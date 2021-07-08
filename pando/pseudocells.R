require(plyr, quietly=T)
require(tidyverse, quietly=T)
require(viridis, quietly=T)
require(gridExtra, quietly=T)
require(dbscan, quietly=T)
require(matrixStats, quietly=T)
require(igraph, quietly=T)
require(ape, quietly=T)
require(reticulate, quietly=T)
require(rdist, quietly=T)
require(ggdendro, quietly=T)
require(tidygraph, quietly=T)
require(ggraph, quietly=T)


invert_list <- function(lov){
    split(rep(names(lov), lengths(lov)), unlist(lov))
}


find_pseudocells <- function(
    object,
    resolution = 0.3,
    k = 10,
    order = 1,
    reduction = 'pca',
    verbose = TRUE,
    return_object = TRUE,
    ...
){

    knn_mat <- FindNeighbors(
        Reductions(object, slot=reduction)@cell.embeddings,
        compute.SNN = FALSE,
        k.param = k,
        verbose = verbose,
        ...
    )$nn

    n_cells <- ncol(knn_mat)
    mask <- runif(n_cells, 0, 1) <= resolution

    g <- graph_from_adjacency_matrix(knn_mat)
    tnn <- ego(g, order=order, nodes=colnames(knn_mat)[mask])
    tnn <- map(tnn, names)
    names(tnn) <- colnames(knn_mat)[mask]
    log_message('Found ', length(tnn), ' territories', verbose=verbose)
    terr <- invert_list(tnn)
    nterr <- mean(map_int(terr, length))
    log_message('On average ', round(nterr, 1), ' territories per cell', verbose=verbose)

    nsing <- 0
    terr_df <- map_dfr(colnames(knn_mat), function(x){
        if (!is.null(terr[[x]])){
            row <- c(x, sample(terr[[x]], 1))
        } else {
            row <- c(x, x)
            nsing <- nsing + 1
        }
        names(row) <- c('cell', 'pseudocell')
        return(row)
    })
    nps <- length(unique(terr_df$pseudocell))
    log_message('Identified ', nps, ' pseudocells',
        ' for ', ncol(knn_mat), ' cells', verbose=verbose)
    if (return_object){
        return(AddMetaData(object, column_to_rownames(terr_df, 'cell')))
    } else {
        return(terr_df)
    }
}
