require(tidyverse, quietly=T)
require(Matrix, quietly=T)
require(sparseMatrixStats, quietly=T)

zscale <- function(x, center=T, scale=T){
    if (center & scale){
        return((x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
    } else if (center & !scale){
        return(x - mean(x, na.rm=TRUE))
    } else if (!center & scale){
        return(x / sd(x, na.rm=TRUE))
    }
}

scale01 <- function(x){
    return((x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))
}

logistic <- function(x){
    return(1/(1+exp(-x)))
}

logoddsratio <- function(x, y=NULL){
    if (!is.null(y)){
        x <- table(x, y)
    }
    lodds <- log((x[1,1]/x[1,2])/(x[2,1]/x[2,2]))
    return(lodds)
}

sparse_cov <- function(X, Y=NULL, XtX=NULL) {
    if (!is(X,"dgCMatrix")) stop("X should be a dgCMatrix")
    if (is.null(Y)) {
        if (is.null(XtX)) {
            XtX <- crossprod(X)
        } else {
            if (ncol(XtX) != ncol(X) | nrow(XtX) != ncol(X)) stop("XtX should have same number of rows and columns as number of columns of X")
        }
        n <- nrow(X)
        cMeans <- colMeans(X)
        covmat <- (as.matrix(XtX) - n*tcrossprod(cMeans))/(n-1)
        sdvec <- sqrt(diag(covmat))
        cormat <- covmat/crossprod(t(sdvec))
        return(list(cov=covmat,cor=cormat))
    } else {
        if (!is(Y,"dgCMatrix")) stop("Y should be a dgCMatrix")
        if (nrow(X) != nrow(Y)) stop("X and Y should have the same number of rows")
        n <- nrow(X)
        cMeansX <- colMeans(X)
        cMeansY <- colMeans(Y)
        covmat <- (as.matrix(crossprod(X,Y)) - n * tcrossprod(cMeansX,cMeansY))/(n-1)
        sdvecX <- sqrt(diag((as.matrix(crossprod(X)) - n*tcrossprod(cMeansX))/(n-1)))
        sdvecY <- sqrt(diag((as.matrix(crossprod(Y)) - n*tcrossprod(cMeansY))/(n-1)))
        cormat <- covmat/outer(sdvecX,sdvecY)
        return(list(cov=covmat,cor=cormat))
    }
}

log_message <- function(..., verbose=T){
    if (verbose){
        message(paste0(...))
    }
}


summary_fun <- list(
    'mean' = sparseMatrixStats::colMeans2,
    'median' = sparseMatrixStats::colMedians,
    'max' = sparseMatrixStats::colMaxs,
    'min' = sparseMatrixStats::colMins,
    'count' = sparseMatrixStats::colCounts,
    'any' = sparseMatrixStats::colAnys,
    'all' = sparseMatrixStats::colAlls,
    'sd' = sparseMatrixStats::colSds,
    'mad' = sparseMatrixStats::colMads
)


aggregate_matrix <- function(
    x,
    groups = NULL,
    fun = 'mean'
){
    if (length(groups) == nrow(x) & 'character'%in%class(fun)){
        if (fun%in%c('count', 'sum')){
            agg_mat <- Matrix.utils::aggregate.Matrix(x=x, groupings=groups, fun=fun)
            return(agg_mat)
        }

        if (fun=='mean'){
            group_counts <- as.numeric(table(groups))
            agg_mat <- Matrix.utils::aggregate.Matrix(x=x, groupings=groups, fun='sum')
            agg_mat <- agg_mat / group_counts
            return(agg_mat)
        }
    }

    if ('character'%in%class(fun)){
        fun <- summary_fun[[fun]]
    }

    if (length(groups) == nrow(x)){
        agg_mat <- sapply(levels(factor(groups)), function(g){
            chunk <- x[which(groups==g), ]
            if (is.null(dim(chunk))){
                return(chunk)
            } else {
                return(fun(chunk))
            }
        })
        agg_mat <- Matrix::Matrix(agg_mat, sparse=T)
    } else if (length(groups) <= 1){
        agg_mat <- fun(x)
        agg_mat <- Matrix::Matrix(agg_mat, sparse=T)
        colnames(agg_mat) <- groups
        rownames(agg_mat) <- colnames(x)
    } else {
        stop('Length of groups must be either nrow(x) or 1.')
    }
    return(Matrix::t(agg_mat))
}


mode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

binarize <- function(expr_mat, thresh=NULL){
    if (!is.null(thresh)){
        expr_mat[expr_mat>thresh] <- 1
    } else {
        expr_mat[expr_mat!=0] <- 1
    }
    return(expr_mat)
}


map_par <- function(x, fun, parallel=FALSE, verbose=TRUE){
    if (!parallel & (verbose==1)){
        return(pbapply::pblapply(X=x, FUN=fun))
    }
    if (!parallel & (verbose!=1)){
        return(base::lapply(X=x, FUN=fun))
    }
    if (parallel){
        outlist <- foreach::foreach(i=1:length(x)) %dopar% {fun(x[[i]])}
        names(outlist) <- names(x)
        return(outlist)
    }
}
