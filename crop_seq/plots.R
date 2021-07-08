#' Plots KO probability distribution for each target gene
#' @rdname plot_ko_probs
#' @export plot_ko_probs
#'
#'
plot_ko_probs <- function(object, ...){
    UseMethod(generic = 'plot_ko_probs', object = object)
}

#' @rdname plot_ko_probs
#' @export
#' @method plot_ko_probs Seurat
#'
#'
plot_ko_probs.Seurat <- function(
    object,
    assay = 'perturb',
    guide_sep = '-'
){

    guide_pattern <- paste0('(.+)', guide_sep, '(\\d+)')

    plot_df <- GetAssayData(object, assay = assay) %>%
        tibble::as_tibble(rownames = 'guide') %>%
        tidyr::gather(cell, p_ko, -guide) %>%
        dplyr::filter(p_ko>0) %>%
        dplyr::mutate(target_gene=str_replace(guide, guide_pattern, '\\1')) %>%
        dplyr::mutate(guide_num=str_replace(guide, guide_pattern, '\\2'))

    p <- ggplot(plot_df, aes(p_ko, guide_num, fill=guide_num)) +
        geom_density_ridges(alpha=0.8) +
        geom_vline(xintercept = 0.5, alpha=0.8) +
        facet_wrap(~target_gene, scales='free') 

    return(p)
}

#' Plots guide assignments on low dim embedding
#' @rdname guide_dim_plot
#' @export guide_dim_plot
#'
#'
guide_dim_plot <- function(object, ...){
    UseMethod(generic = 'guide_dim_plot', object = object)
}

#' @rdname guide_dim_plot
#' @export
#' @method guide_dim_plot Seurat
#'
#'
guide_dim_plot.Seurat <- function(
    object,
    guide_assay = 'guide_assignments',
    reduction = 'umap',
    guide_sep = '-',
    ...
){
    guide_pattern <- paste0('(.+)', guide_sep, '(\\d+)')

    guide_assignments <- t(GetAssayData(object, assay = guide_assay))
    plot_df <- get_guide_meta(
        object,
        guide_assay = guide_assay,
        guide_sep = guide_sep,
        with_reduction = reduction,
    )

    embedding <- Reductions(object, slot = reduction)@cell.embeddings
    object <- AddMetaData(object, as.data.frame(embedding))
    dim1 <- colnames(embedding)[1]
    dim2 <- colnames(embedding)[2]

    guide_num_colors <- c('#C21F5B', '#DE78A4', '#F9E6EE', 'gray')
    names(guide_num_colors) <- c(1,2,3,NA)
    plots <- map(levels(factor(plot_df$target_gene)), function(i){
        p <- ggplot(plot_df, aes(.data[[dim1]], .data[[dim2]], fill=as.character(guide_num))) +
            geom_point(data=dplyr::filter(plot_df, target_gene!=i), color='gray', alpha=0.5) +
            geom_point(data=dplyr::filter(plot_df, target_gene==i), shape=21, stroke=0.1) +
            scale_fill_manual(values=guide_num_colors) +
            theme_void() +
            labs(title=i)
        return(p)
    })

    p <- patchwork::wrap_plots(plots, ...)

    return(p)
}


#' Gets guide df from seurat object
#' @rdname get_guide_meta
#' @export get_guide_meta
#'
#'
get_guide_meta <- function(object, ...){
    UseMethod(generic = 'get_guide_meta', object = object)
}


#' @rdname get_guide_meta
#' @export
#' @method get_guide_meta Seurat
#'
#'
get_guide_meta.Seurat <- function(
    object,
    guide_assay = 'guide_assignments',
    guide_sep = '-',
    with_reduction = 'umap',
    with_meta = TRUE
){

    guide_pattern <- paste0('(.+)', guide_sep, '(\\d+)')

    guide_assignments <- t(GetAssayData(object, assay = guide_assay))
    guide_assign_df <- guide_assignments %>%
        tibble::as_tibble(rownames='cell') %>%
        tidyr::gather(guide, guide_umi_count, -cell) %>%
        dplyr::filter(guide_umi_count!=0) %>%
        dplyr::mutate(target_gene=str_replace(guide, guide_pattern, '\\1')) %>%
        dplyr::mutate(guide_num=str_replace(guide, guide_pattern, '\\2'))

    embedding <- Reductions(object, slot = with_reduction)@cell.embeddings
    object <- AddMetaData(object, as.data.frame(embedding))
    dim1 <- colnames(embedding)[1]
    dim2 <- colnames(embedding)[2]

    guide_assign_df <- object@meta.data %>%
        tibble::as_tibble(rownames='cell') %>%
        dplyr::left_join(guide_assign_df) %>%
        dplyr::mutate(guide_umi_count=ifelse(is.na(guide_umi_count), 0, guide_umi_count)) %>%
        dplyr::mutate(guide_num=ifelse(target_gene=='DUMMY', 1, guide_num)) %>%
        dplyr::group_by(cell) %>%
        dplyr::mutate(n_guides=ifelse(is.na(guide), 0, n())) %>%
        dplyr::ungroup()

    return(guide_assign_df)
}
