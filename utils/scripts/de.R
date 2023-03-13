require(plyr, quietly=T)
require(tidyverse, quietly=T)
require(gridExtra, quietly=T)
require(ggrepel, quietly=T)
require(matrixStats, quietly=T)
require(pbapply, quietly=T)
require(presto, quietly=T)
require(coin, quietly=T)
require(ROCR, quietly=T)
require(pROC, quietly=T)
require(MASS, quietly=T)
require(egg, quietly=T)
require(ggnewscale, quietly=T)

source('~/scripts/single_cell/markers.R')
source('~/scripts/utils/utils.R')
source('~/scripts/utils/colors.R')

#### GLM DE ####
# glm_fams <- list(
#     'gaussian' = gaussian(),
#     'poisson' = poisson()
# )
#
# glm_de <- function(X, Y, family='gaussian', term=1, test='f', adjust_method='holm'){
#     if (test == 'wald'){
#         test_df <- wald_test(X, Y, family=family, term=term)
#     }
#     if (test == 'f'){
#         test_df <- f_test(X, Y, family=family, term=term)
#     }
#     test_df$padj <- p.adjust(test_df$pval, method=adjust_method)
#     return(test_df)
# }


# Optimally, one would interate through terms to get all individual p-values, if I
# write a package, I'll fix that!
# wald_test <- function(X, Y, family='gaussian', term=2){
#     X <- as.matrix(X)
#     tests <- future_apply(Y, 2, function(y){
#         f <- fit_glm(X, y, family=family)
#         vc <- vcov(f) %>% {.[is.na(.)] <- 0; .}
#         cf <- coef(f) %>% {.[is.na(.)] <- 0; .}
#         test <- aod::wald.test(vc, cf, Terms=term)
#         return(list(test=test, coef=coef(f)[term]))
#     })
#     pvals <- map_dbl(tests, function(t){
#         t$test$result$chi2['P']
#     })
#     chi2s <- map_dbl(tests, function(t){
#         t$test$result$chi2['chi2']
#     })
#     coefs <- map_dbl(tests, function(t){
#         t$coef
#     })
#     out_df <- tibble(
#         gene = names(tests),
#         coef = coefs,
#         chi2 = chi2s,
#         pval = pvals
#     )
#     return(out_df)
# }
#
# f_test <- function(X, Y, family='gaussian', term=1){
#     X <- as.matrix(X)
#     tests <- future_apply(Y, 2, function(y){
#         f <- fit_glm(X, y, family=family)
#         name <- names(coef(f))[term]
#         name <- name[!is.na(coef(f)[name])]
#         test <- summary(f)
#         coef <- test$coef[name, 1]
#         std_err <- test$coef[name, 2]
#         tval <- test$coef[name, 3]
#         pval <- test$coef[name, 4]
#         out_df <- tibble(
#             name = name,
#             coef = coef,
#             std_err = std_err,
#             tval = tval,
#             pval = pval
#         )
#         return(out_df)
#     })
#     return(bind_rows(tests, .id='gene'))
# }
#
# fit_glm <- function(X, y, family='gaussian'){
#     if (family %in% c('gaussian', 'poisson')){
#         fit <- glm.fit(X, y, family=glm_fams[[family]])
#         fit$x <- X
#         names(fit$coefficients) <- colnames(X)
#         class(fit) <- c(fit$class, c("glm", "lm"))
#     }
#     if (family == 'negbin'){
#         fit1 <- glm.fit(X, y, family = poisson())
#         theta <- theta.ml(y, mu=fit1$fitted)
#         fit <- glm.fit(
#             X, y,
#             family = MASS::negative.binomial(theta)
#         )
#     }
#     return(fit)
# }

#### UTIL FUNC ####
expression_fraction <- function(expr_mat, margin=2){
    # Percent of non-zero values for each col/row
    if (is.null(dim(expr_mat))){
        return(sum(expr_mat > 0) / length(expr_mat))
    }
    prc_ex <- apply(expr_mat, margin, function(x) {
            return(sum(x > 0) / length(x))
        }
    )
}

percent_expression <- function(expr_mat, groups=NULL, margin=2){
    # Percent of non-zero values for each group
    if (!is.null(groups) & length(levels(factor(groups))) > 1){
        gset <- levels(factor(groups))
        prc_ex <- map(gset, function(x){
            curr_group <- expr_mat[groups == x, ]
            return(expression_fraction(curr_group, margin))
        })
        prc_ex <- do.call(cbind, prc_ex)
        colnames(prc_ex) <- levels(factor(groups))
    } else {
        prc_ex <- expression_fraction(expr_mat)
    }
    return(prc_ex)
}

generalized_fold_change <- function(expr_mat, groups, use_median=F){
    # Computes average expression fold change between groups over quantiles
    if (use_median){
        mean <- median
    }
    gs <- levels(groups)
    fc <- pbapply(expr_mat, 1, function(x){
        q1 <- quantile(x[groups == gs[1]], probs = seq(0.1, 0.9, length.out=9))
        q2 <- quantile(x[groups == gs[2]], probs = seq(0.1, 0.9, length.out=9))
        return(sum(q1 - q2) / 9)
    })
    return(fc)
}


#### Fold change ####
fc <- function(object, ...){
    UseMethod(generic = 'fc', object = object)
}

fc.default <- function(object, groups, use_median=F, na_rm=F){
    # Computes average expression fold change between groups
    if (use_median){
        mean <- median
    }
    gs <- levels(factor(groups))
    fold_change <- map(gs, function(g){
        g1 <- colMeans(object[groups == g, , drop=F], na.rm=na_rm)
        g2 <- colMeans(object[groups != g, , drop=F], na.rm=na_rm)
        fold <- g1 - g2
        out_tbl <- tibble(feature=colnames(object), group=g, fc=fold)
        return(out_tbl)
    })
    return(bind_rows(fold_change))
}

fc.Seurat <- function(object, groups, assay='RNA', slot='data', ...){
    expr_mat <- t(GetAssayData(object, assay=assay, slot=slot))
    groups <- object@meta.data[[groups]]
    return(fc(object=expr_mat, groups=groups, ...))
}

#### Wilcox DE ####
de <- function(object, ...){
    UseMethod(generic = 'de', object = object)
}

de.default <- function(object, groups){
    de_df <- wilcoxauc(t(object), groups)
    de_df <- dplyr::select(de_df, feature, group, 'avg_exp'=avgExpr, 'fc'=logFC, auc, pval, padj, 'prcex_self'=pct_in, 'prcex_other'=pct_out)
    return(as_tibble(de_df))
}

de.Seurat <- function(object, groups, assay='RNA', slot='data'){
    expr_mat <- t(GetAssayData(object, assay=assay, slot=slot))
    groups <- object@meta.data[[groups]]
    return(de(object=expr_mat, groups=groups))
}


nested_de <- function(expr, meta, group1_col, group2_col){
    target_de <- future_lapply(unique(meta[[group1_col]]), function(x){
        group1 <- meta[meta[[group1_col]] == x, ]
        future_lapply(unique(group1[[group2_col]]), function(y){
            print(c(x,y))
            group2_cells <- group1[group1[[group2_col]] == y, ]$cell
            group1_mat <- expr[rownames(expr)%in%group1$cell, ]
            is_g2 <- rownames(group1_mat)%in%group2_cells
            de_df <- de(group1_mat, is_g2) %>%
                dplyr::filter(group==T) %>%
                mutate(group1=x, group=y) %>%
                return()
        }) %>% purrr::reduce(bind_rows)
    }) %>% purrr::reduce(bind_rows)
    return(target_de)
}


grouped_de <- function(expr, meta, group1_col, group2_col){
    target_de <- future_lapply(unique(meta[[group1_col]]), function(x){
        group1 <- meta[meta[[group1_col]] == x, ]
        de_df <- de(expr[group1$cell, ], group1[[group2_col]]) %>%
            mutate(group1=x) %>%
            return()
    }) %>% purrr::reduce(bind_rows)
    return(target_de)
}


#### PLOTS PLOTS PLOTS ####
get_features <- function(object, ...){
    UseMethod(generic = 'get_features', object = object)
}

get_features.default <- function(object, features, gather_output=T){
    data <- as.matrix(object[, colnames(object)%in%features])
    data <- as_tibble(data, rownames="cell")
    if (gather_output){
        data <- tidyr::gather(data, gene, expr, -cell)
    }
    return(data)
}

get_features.Seurat <- function(
    object, features, assay='RNA', slot='data', gather_output=T
){
    data <- GetAssayData(object, assay=assay, slot=slot)
    data <- as.matrix(data[rownames(data)%in%features, ])
    data <- get_features(t(data), features, gather_output=gather_output)
    return(data)
}




point_size <- 0.7
point_alpha <- 0.7
feature_color_scale <- scale_color_gradientn(colours=gyrdpu)
feature_fill_scale <- scale_fill_gradientn(colours=gyrdpu)
feature_theme <- theme(legend.position = 'none')

tsne_plot <- function(x, g){
    ggplot(x, aes(tSNE1, tSNE2, col=expr)) +
        geom_point(size=point_size, alpha=point_alpha) +
        theme_dr() +
        labs(title=g) +
        feature_color_scale +
        feature_theme
}

umap_plot <- function(x, g){
    ggplot(x, aes(UMAP1, UMAP2, col=expr)) +
        geom_point(size=point_size, alpha=point_alpha) +
        theme_dr() +
        labs(title=g) +
        feature_color_scale +
        feature_theme
}

phate_plot <- function(x, g){
    ggplot(x, aes(PHATE1, PHATE2, col=expr)) +
        geom_point(size=point_size, alpha=point_alpha) +
        theme_dr() +
        labs(title=g) +
        feature_color_scale +
        feature_theme
}

fa_plot <- function(x, g){
    ggplot(x, aes(FA1, FA2, col=expr)) +
        geom_point(size=point_size, alpha=point_alpha) +
        theme_dr() +
        labs(title=g) +
        feature_color_scale +
        feature_theme
}

xy_plot <- function(x, g){
    ggplot(x, aes(x, y, fill=expr, alpha=expr)) +
        geom_tile(alpha=point_alpha) +
        theme_dr() +
        labs(title=g) +
        feature_fill_scale +
        feature_theme
}

# feature_plot <- function(expr_mat, meta, features,
#                          plot=umap_plot,
#                          title=NULL,
#                          ncol=NULL,
#                          nrow=NULL,
#                          sort=T,
#                          scale=T){
#     meta$cell <- as.character(meta$cell)
#     expr_features <- get_features(expr_mat, features, scale=scale) %>%
#         dplyr::mutate(cell=as.character(cell)) %>%
#         dplyr::right_join(meta)
#     plots <- map(unique(features), function(g){
#         x <- expr_features %>%
#             dplyr::filter(gene==g)
#         if (sort){
#           x <- dplyr::arrange(x, expr)
#         }
#         plot(x, g)
#     })
#     p_full <- wrap_plots(plots=plots, ncol=ncol, nrow=nrow)
#     return(p_full)
# }

feature_plot <- function(object, colors = gyrdpu(), ...){
    p <- FeaturePlot(object, ...) &
        theme_dr() & theme(axis.line = element_blank()) &
        scale_color_gradientn(colors=colors) & no_legend()
    return(p)
}

dim_plot <- function(object, ...){
    p <- DimPlot(object, ...) &
        theme_void()
    return(p)
}

vulcano_plot <- function(features_df, top=10, top_only=T){

    if (class(top) == 'numeric'){
        top_genes <- features_df %>%
            dplyr::group_by(group) %>%
            dplyr::top_n(top*5, -padj) %>%
            dplyr::group_by(group)
        if (top_only){
            top_genes <- top_genes %>% dplyr::top_n(top, fc)
        } else {
            top_genes <- top_genes %>% dplyr::top_n(top, abs(fc))
        }
        top_genes <- top_genes %>%
            dplyr::mutate(top=TRUE)
    } else {
        top_genes <- features_df %>%
            dplyr::filter(feature%in%top) %>%
            dplyr::mutate(top=TRUE)
    }

    features_df <- features_df %>%
        dplyr::left_join(top_genes) %>%
        dplyr::mutate(x = replace_na(top, FALSE))

    p <- ggplot(features_df, aes(fc, -log10(pval), fill=auc)) +
        geom_point(size=2, shape=21, stroke=0.2) +
        geom_text_repel(
            data = dplyr::filter(features_df, top),
            mapping = aes(fc, -log10(pval), label=feature),
            max.overlaps = 10000

        ) +
        geom_hline(yintercept = -log10(0.01)) +
        scale_fill_gradientn(colors=pals::ocean.curl(100), limits=c(0, 1)) +
        facet_wrap(~group, scales='free') +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size=12, face='bold')
        )
    return(p)
}

vulcauc_plot <- function(features_df, top=10){

    if (class(top) == 'numeric'){
        top_genes <- features_df %>%
            dplyr::group_by(group) %>%
            dplyr::top_n(top*4, -padj) %>%
            dplyr::group_by(group) %>%
            dplyr::top_n(top, fc) %>%
            dplyr::mutate(top=TRUE)
    } else if (!is.null(top)) {
        top_genes <- features_df %>%
            dplyr::filter(feature%in%top) %>%
            dplyr::mutate(top=TRUE)
    }

    if (!is.null(top)){
        features_df <- features_df %>%
            dplyr::left_join(top_genes) %>%
            dplyr::mutate(x = replace_na(top, FALSE))
    } else {
        features_df$top <- FALSE
    }

    p <- ggplot(features_df, aes(auc, -log10(pval), fill=fc)) +
        geom_point(size=2, shape=21, stroke=0.2) +
        geom_text_repel(
            data = dplyr::filter(features_df, top),
            mapping = aes(auc, -log10(pval), label=feature),
            max.overlaps = 10000

        ) +
        geom_hline(yintercept = -log10(0.01)) +
        scale_fill_gradientn(colors=pals::ocean.curl(100), limits=c(0, 1)) +
        facet_wrap(~group, scales='free') +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size=12, face='bold')
        )
    return(p)
}

de_heatmap_plot <- function(features_df, expr_mat, groups,
    top = 10,
    colors = blues_flat,
    group_colors = many,
    newpage = T
){

  if (class(top) == 'numeric'){
    top_features <- features_df %>%
      filter(padj<1e-5, auc>0.5) %>%
      group_by(group) %>%
      top_n(top + 20, fc) %>%
      ungroup() %>%
      group_by(feature) %>%
      dplyr::mutate(max_auc = max(auc),
             max_fc = max(fc)) %>%
      filter(fc==max_fc) %>%
      ungroup() %>%
      group_by(group) %>%
      top_n(top, fc) %>%
      arrange(group, desc(fc)) %>%
      {unique(.$feature)}
  } else {
    top_features <- unique(top)
  }

  top_marker_expr <- get_features(expr_mat, top_features) %>%
    mutate(groups=factor(groups[match(.$cell, rownames(expr_mat))], levels=levels(factor(groups)))) %>%
    group_by(groups) %>%
    mutate(mean_expr=mean(expr)) %>%
    ungroup() %>%
    arrange(groups, mean_expr) %>%
    mutate(cell=factor(cell, levels=unique(cell)), feature=factor(feature, levels=top_features))

  p1 <- ggplot(top_marker_expr, aes(cell, feature, fill=expr)) +
    scale_fill_gradientn(colors=colors) +
    geom_raster() +
    scale_y_discrete(expand=expand_scale(add = c(0,0))) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none'
    )
  p2 <- ggplot(top_marker_expr, aes(cell, 1, fill=groups)) +
    geom_raster() +
    scale_fill_manual(values=group_colors) +
    scale_y_discrete(expand=expand_scale(add = c(0,0)), labels=function(breaks) {rep_along(breaks, " ")}) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        legend.position="bottom"
    )
  egg::ggarrange(p1, p2, padding=0, heights = c(40, 2), newpage = newpage)
}



heatmap_plot <- function(
    expr_mat,
    groups =  NULL,
    genes = NULL,
    colors = blues_flat,
    group_colors = many,
    newpage = T
){

    if (is.null(genes)){
        top_features <- rownames(expr_mat)
    } else {
        top_features <- unique(genes)
    }

    if (is.null(groups)){
        top_features <- colnames(expr_mat)
    } else {
        top_features <- levels(factor(groups))
    }

    top_marker_expr <- get_features(expr_mat, genes) %>%
        mutate(groups=factor(groups[match(.$cell, rownames(expr_mat))], levels=levels(factor(groups)))) %>%
        group_by(groups) %>%
        mutate(mean_expr=mean(expr)) %>%
        ungroup() %>%
        arrange(groups, mean_expr) %>%
        mutate(cell=factor(cell, levels=unique(cell)),
               gene=factor(gene, levels=genes))

    p1 <- ggplot(top_marker_expr, aes(cell, gene, fill=expr)) +
        scale_fill_gradientn(colors=colors) +
        geom_raster() +
        scale_y_discrete(expand=expand_scale(add = c(0,0))) +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank(),
            legend.position = 'none'
        )
    p2 <- ggplot(top_marker_expr, aes(cell, 1, fill=groups)) +
        geom_raster() +
        scale_fill_manual(values=group_colors) +
        scale_y_discrete(expand=expand_scale(add = c(0,0)), labels=function(breaks) {rep_along(breaks, " ")}) +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            panel.border = element_blank(),
            legend.position="bottom"
        )
    egg::ggarrange(p1, p2, padding=0, heights = c(40, 2), newpage = newpage)

}

#### Enrichment tests ####
# fisher_test_single <- function(x, y, ...){
#     x_out <- c()
#     y_out <- c()
#     orats <- c()
#     pvals <- c()
#     for (x_ in unique(x)){
#         for (y_ in unique(y)){
#             if (length(unique(x)) > 1 & length(unique(y)) > 1){
#                 ftest <- fisher.test(table(x==x_, y==y_), ...)
#                 orats <- c(orats, ftest$estimate)
#                 pvals <- c(pvals, ftest$p.value)
#                 x_out <- c(x_out, x_)
#                 y_out <- c(y_out, y_)
#             } else {
#                 orats <- c(orats, Inf)
#                 pvals <- c(pvals, 1)
#                 x_out <- c(x_out, x_)
#                 y_out <- c(y_out, y_)
#             }
#         }
#     }
#     out_tbl <- tibble(
#         x = x_out,
#         y = y_out,
#         odds_ratio = orats,
#         pval = pvals
#     )
#     return(out_tbl)
# }

# cmh_test_single <- function(x, y, blocks=NULL, ...){
#     x_out <- c()
#     y_out <- c()
#     stats <- c()
#     pvals <- c()
#     for (x_ in unique(x)){
#         for (y_ in unique(y)){
#             x_fac = factor(x==x_)
#             y_fac = factor(y==y_)
#             if (is.null(blocks)){
#                 ftest <- coin::cmh_test(x_fac ~ y_fac, ...)
#             } else {
#                 ftest <- coin::cmh_test(x_fac ~ y_fac | factor(blocks), ...)
#             }
#             stats <- c(stats, statistic(ftest))
#             pvals <- c(pvals, pvalue(ftest))
#             x_out <- c(x_out, x_)
#             y_out <- c(y_out, y_)
#         }
#     }
#     out_tbl <- tibble(
#         x = x_out,
#         y = y_out,
#         statistic = stats,
#         pval = pvals
#     )
#     return(out_tbl)
# }

# grouped_fisher_test <- function(x, y, groups, ...){
#     ftest <- map(unique(groups), function(g){
#         select <- groups == g
#         fisher_test_single(x[select], y[select], ...)
#     }) %>% {names(.) <- unique(groups); .} %>%
#     bind_rows(.id = 'group') %>%
#     mutate(log_odds_ratio = ifelse(is.infinite(log(odds_ratio)), 0, log(odds_ratio)))
# }

# fisher_test <- function(x, y, groups=NULL, ...){
#     if (is.null(groups)){
#         return(fisher_test_single(x, y, ...))
#     } else {
#         return(grouped_fisher_test(x, y, groups))
#     }
# }

# proportion_test <- function(x, y, groups=NULL, blocks=NULL, method='fisher'){
#     if (method == 'fisher'){
#         return(fisher_test(x, y, groups=groups))
#     } else if (method == 'cmh'){
#         return(cmh_test_single(x, y, blocks=blocks))
#     }
# }
