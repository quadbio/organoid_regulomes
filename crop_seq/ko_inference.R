#' Adds guide assay to seurat object
#' @rdname add_guide_assay
#' @export add_guide_assay
#'
#'
add_guide_assay <- function(object, ...){
    UseMethod(generic = 'add_guide_assay', object = object)
}

#' @rdname add_guide_assay
#' @export
#' @method add_guide_assay Seurat
#'
#'
add_guide_assay.Seurat <- function(
    object,
    guide_assignments,
    assay_name = 'guide_assignments'
){
    cells_use <- intersect(colnames(object), rownames(guide_assignments))
    guide_mat <- guide_assignments[cells_use, ]

    # Make probs the right shape
    guide_assay <- matrix(nrow = ncol(object), ncol = ncol(guide_mat))
    guide_assay[is.na(guide_assay)] <- 0
    rownames(guide_assay) <- colnames(object)
    colnames(guide_assay) <- colnames(guide_mat)
    guide_assay <- Matrix::Matrix(guide_assay, sparse=TRUE)
    guide_assay[rownames(guide_mat), colnames(guide_mat)] <- guide_mat

    object[[assay_name]] <- CreateAssayObject(
        data = t(guide_assay)
    )

    return(object)
}


#' Infers KO probs from guide assignments and transcriptome
#' @rdname infer_ko_probs
#' @export infer_ko_probs
#'
#'
infer_ko_probs <- function(object, ...){
    UseMethod(generic = 'infer_ko_probs', object = object)
}


#' @rdname infer_ko_probs
#' @export
#' @method infer_ko_probs Seurat
#'
#'
infer_ko_probs.Seurat <- function(
    object,
    covariates = NULL,
    guide_sep = '-',
    assay = 'RNA',
    slot = 'data',
    alpha = 0.5,
    assay_name = 'perturb',
    guide_assay = 'guide_assignments',
    save_coefs = FALSE,
    genes_use = NULL,
    verbose = TRUE
){
    guide_assignments <- t(GetAssayData(object, assay = guide_assay))
    guide_pattern <- paste0('(.+)', guide_sep, '(\\d+)')
    target_genes <- str_replace(colnames(guide_assignments), guide_pattern, '\\1')
    colnames(guide_assignments) <- str_replace_all(
        colnames(guide_assignments), '-', '_')

    expr_mat <- t(GetAssayData(object, assay = assay, slot = slot))
    cells_use <- intersect(rownames(expr_mat), rownames(guide_assignments))

    X <- binarize(guide_assignments[cells_use, ])
    C <- object@meta.data[cells_use, covariates]

    if (is.null(genes_use)){
        genes_use <- VariableFeatures(object)
    }
    Y <- expr_mat[cells_use, genes_use]

    ko_inf <- infer_ko_probs(
        X, Y, C,
        alpha = alpha,
        target_genes = target_genes,
        verbose = verbose
    )

    # Make probs the right shape
    ko_probs <- matrix(nrow = ncol(object), ncol = ncol(guide_assignments))
    ko_probs[is.na(ko_probs)] <- 0
    rownames(ko_probs) <- colnames(object)
    colnames(ko_probs) <- colnames(guide_assignments)
    ko_probs <- Matrix::Matrix(ko_probs, sparse=TRUE)
    ko_probs[rownames(ko_inf$ko_probabilities), colnames(ko_inf$ko_probabilities)] <- ko_inf$ko_probabilities

    # Save to seurat object: data -> infered probs, counts -> raw assignments
    object[[assay_name]] <- CreateAssayObject(
        data = t(ko_probs)
    )

    object[['guide_assignments']] <- CreateAssayObject(
        data = binarize(t(ko_probs))
    )

    if (save_coefs){
        slot_name = paste0(assay_name, '_coefs')
        Misc(object = object, slot = slot_name) <- ko_inf$coefficients
    }

    return(object)
}


#' @rdname infer_ko_probs
#' @export
#' @method infer_ko_probs default
#'
#'
infer_ko_probs <- function(
    X, Y,
    C = NULL,
    alpha = 0.5,
    family = 'gaussian',
    target_genes = NULL,
    return_models = F,
    parallel = FALSE,
    verbose = T,
    ...
){

    log_message('Preparing data.', verbose=verbose)
    n_genes <- length(unique(target_genes))
    guide_formula <- reformulate(termlabels = c(colnames(C), colnames(X)))

    if (!is.null(target_genes)){
        X_guide <- X
        X <- summarize_target_genes(X_guide, target_genes)
        gene_formula <- reformulate(termlabels = c(colnames(C), colnames(X)))
    } else {
        X_guide <- X
        gene_formula <- guide_formula
    }

    # Leave out first column of the model matrix,
    # as glmnet will add the intercept itself
    log_message('Preparing model input.', verbose=verbose)
    XC <- model.matrix(gene_formula, cbind(C, X))[, -1]
    XC_guide <- model.matrix(guide_formula, cbind(C, X_guide))[, -1]
    Y_sigmas <- sparseMatrixStats::rowSds(Y)

    log_message('Fitting model.', verbose=verbose)
    fits <- fit_glmnet(
        XC, Y, alpha=alpha, family=family, parallel=parallel, verbose=verbose)
    Y <- Y[, names(fits)]

    log_message('Getting coefficients.', verbose=verbose)
    # Cbind later
    coefs <- map_par(fits, get_coef, parallel=parallel, verbose=verbose)
    names(coefs) <- names(fits)
    coefs_guide <- lapply(coefs, function(c){
        c_idx <- c(rownames(c)[1:(nrow(c)-n_genes)], target_genes)
        c_new <- Matrix(c[c_idx, ], sparse=T)
        rownames(c_new) <- c_idx
        rownames(c_new)[(nrow(c_new)-ncol(X_guide)+1):nrow(c_new)] <- colnames(X_guide)
        return(c_new)
    })

    log_message('Calculating error.', verbose=verbose)
    Y_hat <- do.call(
        cbind,
        map_par(coefs_guide, function(m){
            y_pred <- coef_predict(m, XC_guide)
            return(y_pred)
        }, parallel=parallel, verbose=verbose)
    )
    XC_err <- (Y - Y_hat)^2

    log_message('Calculating KO probability.', verbose=verbose)
    ko_prob <- map_par(seq(ncol(X_guide)), function(j){
        X_0 <- X_guide
        X_0[, j] <- 0
        XC_0 <- model.matrix(guide_formula, cbind(C, X_0))[, -1]
        Y_pred <- do.call(
            cbind,
            lapply(coefs_guide, function(m){
                y_pred <- coef_predict(m, XC_0)
                return(y_pred)
            })
        )
        err_diff <- rowSums((Y - Y_pred)^2 - XC_err)
        p_ko <- logistic(err_diff / 2*Y_sigmas**2)
        return(p_ko)
    }, parallel=parallel, verbose=verbose)
    ko_prob <- do.call(cbind, ko_prob)
    ko_prob <- ko_prob * X_guide
    colnames(ko_prob) <- colnames(X_guide)
    ko_prob <- Matrix(ko_prob, sparse=T)

    if (!return_models){
        fits <- NULL
    }

    coef_mat <- do.call(cbind, coefs)
    colnames(coef_mat) <- names(coefs)

    out_list <- list(
        models = fits,
        coefficients = Matrix(coef_mat, sparse=T),
        squared_error = Matrix(XC_err, sparse=T),
        ko_probabilities = ko_prob
    )

    return(out_list)
}


fit_glmnet <- function(X, Y, alpha=0.5, family='gaussian', parallel=F, verbose=T){
    fits <- map_par(seq(ncol(Y)), function(i){
        try(cv.glmnet(
            X, Y[,i],
            alpha = alpha,
            family = family
        ), silent=TRUE)
    }, parallel=parallel, verbose=verbose)
    names(fits) <- colnames(Y)
    fits <- fits[map_chr(fits, class) == 'cv.glmnet']
    return(fits)
}

get_coef <- function(fit){
    if (class(fit) == 'cv.glmnet'){
        return(coef(fit, s='lambda.min'))
    } else if (class(fit) == 'cv.glmreg'){
        coefs <- coef(fit)
        coef_mat <- Matrix(coefs, sparse=TRUE, ncol=1)
        rownames(coef_mat) <- names(coefs)
        return(coef_mat)
    }
}

add_intercept <- function(x){
    incpt <- rep(1,nrow(x))
    x <- as.data.frame(cbind(incpt,x))
    colnames(x)[1] <- '(Intercept)'
    return(x)
}

coef_predict <- function(coefs, X){
    t(tcrossprod(t(coefs), add_intercept(X)))
}

summarize_target_genes <- function(x, target_genes){
    t(binarize(aggregate_matrix(
        t(x),
        groups = target_genes,
        fun = colMeans
    )))
}
