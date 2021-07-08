suppressPackageStartupMessages(library(glmnetUtils, quietly=T, verbose=F))
suppressPackageStartupMessages(library(tibble, quietly=T, verbose=F))


fit_model <- function(
    formula, data,
    method = 'glm',
    family = gaussian,
    alpha = 1,
    ...
){
    result <- switch(
        method,
        'glm' = fit_glm(formula, data, family=family, ...),
        'glmnet' = fit_glmnet(formula, data, family=family, alpha=alpha, ...),
        'cv.glmnet' = fit_cvglmnet(formula, data, family=family, alpha=alpha, ...)
    )
    return(result)
}

fit_glm <- function(formula, data, family=gaussian, ...){
    fit <- glm(formula, data=data, family=family, ...)
    s <- summary(fit)
    gof <- tibble(
        dsq = with(s, 1 - deviance/null.deviance)
    )
    coefs <- as_tibble(s$coefficients, rownames='term')
    colnames(coefs) <- c('term', 'estimate', 'std_err', 'statistic', 'pval')
    return(list(gof=gof, coefs=coefs))
}

fit_glmnet <- function(
    formula, data,
    family = gaussian,
    alpha = 1,
    nlambda = 20,
    ...
){
    fit <- glmnetUtils::glmnet(
        formula,
        data = data,
        family = family,
        alpha = alpha,
        nlambda = nlambda,
        ...
    )
    class(fit) <- 'glmnet'
    which_max <- which(fit$dev.ratio > max(fit$dev.ratio) * 0.95)[1]
    lambda_choose <- fit$lambda[which_max]
    gof <- tibble(
        lambda = lambda_choose,
        dsq = fit$dev.ratio[which_max],
        alpha = alpha
    )
    coefs <- as_tibble(as.matrix(coef(fit, s=lambda_choose)), rownames='term')
    colnames(coefs) <- c('term', 'estimate')
    return(list(gof=gof, coefs=coefs))
}

fit_cvglmnet <- function(
    formula, data,
    family = gaussian,
    alpha = 0.5,
    nlambda = 20,
    nfolds = 5,
    ...
){
    fit <- glmnetUtils::cv.glmnet(
        formula,
        data = data,
        family = family,
        alpha = alpha,
        nlambda = nlambda,
        nfolds = nfolds,
        ...
    )
    class(fit) <- 'cv.glmnet'
    which_max <- fit$index['1se', ]
    gof <- tibble(
        lambda = fit$lambda.1se,
        dsq = fit$glmnet.fit$dev.ratio[which_max],
        alpha = alpha
    )
    coefs <- as_tibble(as.matrix(coef(fit)), rownames='term')
    colnames(coefs) <- c('term', 'estimate')
    return(list(gof=gof, coefs=coefs))
}
