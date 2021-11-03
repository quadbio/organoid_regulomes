require(tidyverse, quietly=T)
require(Matrix, quietly=T)
require(sparseMatrixStats, quietly=T)
require(glmnet, quietly=T)
require(Seurat, quietly=T)
require(lmtest, quietly=T)


#' @rdname glm_de
#' @export glm_de
#' @method glm_de Seurat
glm_de.Seurat <- function(
	object,
	test_var,
	covariates = NULL,
	assay = 'RNA',
	slot = 'data',
	family = 'gaussian',
	test = 't',
	genes_use = NULL,
	adjust_method = 'holm'
){

	expr_mat <- t(GetAssayData(object, assay = assay, slot = slot))
	X <- object@meta.data[, test_var]
	C <- object@meta.data[, covariates, drop=FALSE]

	model_formula <- reformulate(termlabels = c(covariates, 'X'))
	XC <- cbind(C, X)
	XC_model <- model.matrix(model_formula, XC)

	if (!is.null(genes_use)){
		expr_mat <- expr_mat[, genes_use]
	}

	de_df <- glm_de(
		XC_model, expr_mat,
		term = ncol(XC_model),
		family = family,
		test = test,
		adjust_method = adjust_method
	)

	return(de_df)
}


#' Infers KO probs from guide assignments and transcriptome
#' @rdname get_perturbation_effects
#' @export get_perturbation_effects
get_perturbation_effects <- function(object, ...){
    UseMethod(generic = 'get_perturbation_effects', object = object)
}

#' @rdname get_perturbation_effects
#' @export
#' @method get_perturbation_effects Seurat
get_perturbation_effects.Seurat <- function(
    object,
	assay = 'RNA',
	slot = 'data',
    perturb_assay = 'perturb',
	covariates = NULL,
    guide_sep = '-',
	genes_use = NULL,
	family = 'gaussian',
	adjust_method = 'holm',
	bin_thresh = NULL,
	verbose = T,
	parallel = T
){
	guide_pattern <- paste0('(.+)', guide_sep, '(\\d+)')
	ko_prob <- GetAssayData(object, assay = perturb_assay, slot = 'data')
	expr_mat <- t(GetAssayData(object, assay = assay, slot = slot))

	target_genes <- str_replace(rownames(ko_prob), guide_pattern, '\\1')
	X_new <- t(aggregate_matrix(as.matrix(ko_prob), groups = target_genes, fun = colMaxs))
	rownames(X_new) <- colnames(ko_prob)

	if (!is.null(bin_thresh)){
		X_new <- binarize(X_new, thresh=bin_thresh)
	}

	C <- object@meta.data[, covariates]
	XC_new <- cbind(C, X_new)
	model_formula <- reformulate(termlabels = c(covariates, colnames(X_new)))
	XC_model_new <- model.matrix(model_formula, XC_new)
	test_terms <- which(colnames(XC_model_new) %in% colnames(X_new))

	if (is.null(genes_use)){
		genes_use <- VariableFeatures(object, assay = assay)
	}
	Y <- expr_mat[, genes_use]

	refit_de <- glm_de(
		XC_model_new, Y,
		term = test_terms,
		family = family,
		test = 't',
		adjust_method = adjust_method,
		verbose = verbose,
		parallel = parallel
	)

	return(refit_de)
}

#### GLM DE ####
glm_fams <- list(
    'gaussian' = gaussian(),
    'poisson' = poisson()
)


#' Infers KO probs from guide assignments and transcriptome
#' @rdname glm_de
#' @export glm_de
glm_de <- function(object, ...){
    UseMethod(generic = 'glm_de', object = object)
}


#' @rdname glm_de
#' @export glm_de
#' @method glm_de default
glm_de.default <- function(
	object, Y,
	family = 'gaussian',
	term = 1,
	test = 't',
	adjust_method = 'holm',
	verbose = T,
	parallel = T
){
	X <- object
    if (test == 'wald'){
        test_df <- wald_test(X, Y, family=family, term=term)
    }
    if (test == 't'){
        test_df <- t_test(
			X, Y, family=family, term=term, verbose = verbose, parallel = parallel)
    }
    test_df$padj <- p.adjust(test_df$pval, method=adjust_method)
    return(test_df)
}

t_test <- function(X, Y, family='gaussian', term=1, verbose = T, parallel = T){
    X <- as.matrix(X)
    tests <- map_par(1:ncol(Y), function(i){
		out_df <- t_test_(X, Y[,i], family, term)
		if ('tbl' %in% class(out_df)){
			return(out_df)
		}
    }, verbose=verbose, parallel=parallel)
	names(tests) <- colnames(Y)
    return(bind_rows(tests, .id='gene'))
}

t_test_ <- function(X, y, family, term){
	f <- fit_glm(X, y, family=family)
	name <- names(coef(f))[term]
	name <- name[!is.na(coef(f)[name])]
	test <- summary(f)
	coef <- test$coef[name, 1]
	std_err <- test$coef[name, 2]
	tval <- test$coef[name, 3]
	pval <- test$coef[name, 4]
	out_df <- tibble(
		name = name,
		coef = coef,
		std_err = std_err,
		tval = tval,
		pval = pval
	)
	return(out_df)
}

fit_glm <- function(X, y, family='gaussian'){
    if (family %in% c('gaussian', 'poisson')){
        fit <- glm.fit(X, y, family=glm_fams[[family]])
        fit$x <- X
        names(fit$coefficients) <- colnames(X)
        class(fit) <- c(fit$class, c("glm", "lm"))
    }
    if (family == 'negbin'){
        fit1 <- glm.fit(X, y, family = poisson())
        theta <- theta.ml(y, mu=fit1$fitted)
        fit <- glm.fit(
            X, y,
            family = MASS::negative.binomial(theta)
        )
    }
    return(fit)
}


# Optimally, one would interate through terms to get all individual p-values, if I
# write a package, I'll fix that!
wald_test <- function(X, Y, family='gaussian', term=2){
    X <- as.matrix(X)
    tests <- apply(Y, 2, function(y){
        f <- fit_glm(X, y, family=family)
        vc <- vcov(f) %>% {.[is.na(.)] <- 0; .}
        cf <- coef(f) %>% {.[is.na(.)] <- 0; .}
        test <- try(aod::wald.test(vc, cf, Terms=term), silent=T)
		if (any(class(test)=='try-error')){
			return()
		}
        return(list(test=test, coef=coef(f)[term]))
    })
	tests <- tests[!map_lgl(tests, is.null)]
    pvals <- map_dbl(tests, function(t){
        t$test$result$chi2['P']
    })
    chi2s <- map_dbl(tests, function(t){
        t$test$result$chi2['chi2']
    })
    coefs <- map_dbl(tests, function(t){
        t$coef
    })
    out_df <- tibble(
        gene = names(tests),
        coef = coefs,
        chi2 = chi2s,
        pval = pvals
    )
    return(out_df)
}


#' Performs diff expression based on a likelihood ratio test
#' @rdname lr_de
#' @export lr_de
lr_de <- function(object, ...){
  UseMethod(generic = 'lr_de', object = object)
}


#' @rdname lr_de
#' @export lr_de
#' @method lr_de Seurat
lr_de.Seurat <- function(
    object,
    test_var,
    covariates = NULL,
    assay = 'RNA',
    slot = 'data',
    family = 'gaussian',
    adjust_method = 'holm',
    features_use = NULL,
    test_use = VariableFeatures(object, assay=assay)
){
    expr_mat <- t(GetAssayData(object, assay = assay, slot = slot))
    X <- object@meta.data[, test_var]
    C <- object@meta.data[, covariates, drop=FALSE]
    XC <- cbind(C, GROUP=X)

    if (!is.null(features_use)){
      expr_mat <- expr_mat[, features_use]
    }

    de_df <- lr_de(
      XC, expr_mat,
      term = ncol(XC),
      family = family,
      adjust_method = adjust_method
    )

    return(de_df)
}

#' @rdname lr_de
#' @export lr_de
#' @method lr_de default
lr_de.default <- function(
    object, Y,
    family = 'gaussian',
    term = 1,
    adjust_method = 'holm'
){
    X <- object
    test_list <- foreach (n=1:ncol(Y)) %dopar% {
        feature <- colnames(Y)[n]
        y <- as.numeric(Y[, n])

        model_data <- cbind(FEATURE=y, X)

        frml1 <- as.formula(object = paste(
            'FEATURE ~ ', paste(colnames(X), collapse = '+')
        ))
        frml0 <- as.formula(object = paste(
            'FEATURE ~ ', paste(colnames(X[, -term, drop=F]), collapse = '+')
        ))

        fit1 <- glm(frml1, data=model_data, family='binomial')
        fit0 <- glm(frml0, data=model_data, family='binomial')

        test <- lmtest::lrtest(fit1, fit0)
        pval <- test$Pr[2]
        coef <- fit1$coefficients[[term+1]]

        return(tibble(feature=feature, coef=coef, pval=pval))
    }
    return(bind_rows(test_list))
}
