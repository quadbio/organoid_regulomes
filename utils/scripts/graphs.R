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

# fnc
knn_bruteforce <- function(dst, k, margin=1, is_prox=F){
    nn_mat <- apply(as.matrix(dst), margin, function(x){
        nn <- which(x <= sort(x, decreasing=is_prox)[k])
        nn <- nn[1:(k)]
        names(nn) <- NULL
        return(nn)
    })
    knn_tbl <- gather(as_tibble(t(nn_mat), rownames="source"),
            key=id, value=idx, -source) %>%
        select(source, idx) %>%
        mutate(target=rownames(as.matrix(dst))[.$idx]) %>%
        select(source, target)
    return(knn_tbl)
}

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

knn_drc <- function(dst, k, grp, drc){
    mtx <- as.matrix(dst)
    names(grp) <- rownames(mtx)
    knn_tbl <- map_dfr(rownames(mtx), function(x){
        drc_idx <- which(drc==grp[x])
        grps_allow <- drc[max(1, drc_idx-1) : min(length(drc), drc_idx+1)]
        nn_allow <- mtx[x,][names(grp)[grp%in%grps_allow]]
        nn <- names(sort(nn_allow)[2:(k+1)])
        dst <- sort(nn_allow)[2:(k+1)]
        data_frame(source=rep(x, k), target=nn, value=dst)
    })
    return(knn_tbl)
}

nn2mat <- function(nn, offset1, offset2, dims) {
    k <- ncol(x = nn)
    j <- as.numeric(x = t(x = nn)) + offset2
    i <- ((1:length(x = j)) - 1) %/% k + 1 + offset1
    nn.mat <- sparseMatrix(i = i, j = j, x = 1, dims = dims)
    return(nn.mat)
}


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


as_node_df <- function(adjacency_df, layout=nicely()){
    adj_graph <- graph_from_data_frame(adjacency_df, directed=T)
    lay <- layout_(adj_graph, layout=nicely())
    adjacency_df <- adjacency_df[, 1:2]
    colnames(adjacency_df) <- c('node1', 'node2')

    node_df <- as_tibble(lay) %>%
        mutate(node2=vertex_attr(adj_graph)$name) %>%
        {colnames(.)[1:2] <- c('x', 'y'); .} %>%
        left_join(adjacency_df) %>%
        mutate(node1=ifelse(is.na(node1), node2, node1))

    from_pos <- filter(node_df, node2%in%adjacency_df[[1]]) %>%
        {.[, 1:3]} %>%
        {colnames(.)[1:3] <- c('x_from', 'y_from', 'node1'); .} %>%
        distinct(x_from, y_from, node1)

    to_pos <- filter(node_df, node2%in%adjacency_df[[2]]) %>%
        {colnames(.)[1:3] <- c('x_to', 'y_to', 'node2'); .} %>%
        distinct(x_to, y_to, node2)

    edge_df <- inner_join(node_df, from_pos) %>%
        left_join(to_pos) %>%
        mutate(
            x_to=ifelse(is.na(x_to), x_from, x_to),
            y_to=ifelse(is.na(y_to), y_from, y_to)) %>%
        select(node1, node2, x, y, x_from, y_from, x_to, y_to)

    return(edge_df)
}


as_dendro_df <- function(
    x,
    dist_method = 'euclidean',
    clust_method = 'ward.D2'
){
    if (class(x) != 'dist'){
        x <- dist(x, method=dist_method)
    } else {
        x <- as.dist(x)
    }
    x <- hclust(x, method=clust_method)
    x <- as.dendrogram(x)
    x <- segment(dendro_data(x))
    return(x)
}
