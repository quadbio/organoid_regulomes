
suppressMessages(source('pando/atac.R'))
suppressMessages(source('pando/models.R'))
suppressMessages(source('utils/utils.R'))

suppressMessages(library(sparseMatrixStats, quietly=T, verbose=F))
suppressMessages(library(Matrix.utils, quietly=T, verbose=F))
suppressPackageStartupMessages(library(optparse, quietly=T, verbose=F))
suppressPackageStartupMessages(library(pbapply, quietly=T, verbose=F))
suppressPackageStartupMessages(library(glmnetUtils, quietly=T, verbose=F))

### CONST
GENE_ANNOT <- 'EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds'
MOTIF2TF <- 'motif2tf_early_all.tsv'

### OPTS
parser <- OptionParser()

parser <- add_option(parser,
    c('-i', '--input'),
    dest='input',
    default='./',
    help='Seurat object as rds file.')

parser <- add_option(parser,
    c('-o', '--out-dir'),
    dest='out',
    default='./',
    help='Output directory.')

parser <- add_option(parser,
    c('-d', '--upstream-dist'),
    dest='upstream_dist',
    default='10000',
    type='integer',
    help='Distance upstream of TSS to consider.')

parser <- add_option(parser,
    c('-t', '--interaction'),
    dest='inter_term',
    default=':',
    type='character',
    help='Interaction term for peak and transcription factor.')

parser <- add_option(parser,
    c('--tf-cor'),
    dest='tf_cor',
    default='0.1',
    type='numeric',
    help='Threshold for TF - target gene correlation.')

parser <- add_option(parser,
    c('--peak-cor'),
    dest='peak_cor',
    default='0',
    type='numeric',
    help='Threshold for peak - target gene correlation.')

parser <- add_option(parser,
    c('-s', '--peak-subset'),
    dest='peak_subset',
    type='character',
    help='Selection criterion for peaks.')

parser <- add_option(parser,
    c('-m', '--method'),
    dest='method',
    type='character',
    default='glm',
    help='Method to fit glm.')

parser <- add_option(parser,
    c('-a', '--alpha'),
    dest='alpha',
    type='numeric',
    default='0.5',
    help='Alpha for glm regularization.')

parser <- add_option(parser,
    c('-c', '--n-cores'),
    dest='cores',
    default='1',
    type='integer',
    help='Number of cpus.')

parser <- add_option(parser,
    c('--verbose'),
    dest='verbose',
    default='F',
    action='store_true',
    help='Verbosity.')


### MAIN
args <- parse_args(parser)
if (is.null(args$verbose)){
    verbose <- FALSE
} else {
    verbose <- args$verbose
}
if (args$cores>1){
    options(future.globals.maxSize = 100000 * 1024^2)
    options(future.rng.onMisuse="ignore")
    plan('multicore', workers=args$cores)
    lapply_use <- future_lapply
    print_genes <- verbose
} else {
    lapply_use <- pblapply
    print_genes <- FALSE
}

# Read input
log_message('Reading input.', verbose=verbose)
srt <- read_rds(args$input)
gene_annot <- read_rds(GENE_ANNOT)
motif2tf <- suppressMessages(read_tsv(MOTIF2TF)) %>%
    dplyr::select(motif, tf) %>%
    distinct() %>%
    mutate(val=1) %>%
    pivot_wider(names_from = 'tf', values_from=val, values_fill=0) %>%
    as_matrix() %>% Matrix::Matrix(sparse=T)

dir.create(args$out, recursive=TRUE)
setwd(args$out)

# Select genes
log_message('Selecting genes.', verbose=verbose)
srt@active.assay <- 'RNA'
gene_data <- GetAssayData(srt, assay='RNA', slot='data')
genes_use <- srt@assays$RNA@misc$features[['var_features']] %>%
    intersect(rownames(gene_data)) %>% intersect(gene_annot$gene_name)
tfs_use <- srt@assays$RNA@misc$features[['var_tfs']] %>%
    intersect(rownames(gene_data)) %>% intersect(gene_annot$gene_name)
gene_annot_use <- gene_annot[gene_annot$gene_name%in%genes_use, ]


# Select peaks at gene
log_message('Selecting peaks at genes.', verbose=verbose)
srt@active.assay <- 'peaks'
peak_data <- srt@assays$peaks@misc$agg$louvain_100
cluster_groups <- srt$peaks_snn_res.100

if (!is.null(args$peak_subset)){
    sub_peaks <- srt@assays$peaks@misc$subsets[[args$peak_subset]]
    peak_ranges <- sub_peaks$ranges
    peak_data <- peak_data[, sub_peaks$peaks]
    colnames(peak_data) <- rownames(sub_peaks$motifs@data)
    peaks2motif <- sub_peaks$motifs@data
} else {
    peak_ranges <- StringToGRanges(rownames(peak_data))
    peaks2motif <- Motifs(srt, assay='peaks')@data
}

peaks_near_gene <- find_peaks_near_genes(
    peak_ranges, gene_annot_use, distance=args$upstream_dist, only_tss=F)
peaks2gene <- aggregate.Matrix(t(peaks_near_gene), groupings=colnames(peaks_near_gene))
peaks_at_gene <- as.logical(sparseMatrixStats::colMaxs(peaks2gene))


# Select motifs in peaks
log_message('Selecting transcription factors with motifs in peaks.', verbose=verbose)
peaks_with_motif <- as.logical(sparseMatrixStats::rowMaxs(peaks2motif))
good_peaks <- colMaxs(peak_data) > 0.1
peaks_use <- peaks_at_gene & peaks_with_motif & good_peaks
peaks2gene <- peaks2gene[, peaks_use]
peaks2motif <- peaks2motif[peaks_use, ]
peak_data <- peak_data[, peaks_use]

# Prepare model input
gene_data <- t(gene_data)
cluster_groups <- cluster_groups[rownames(gene_data)]
tfs_use <- intersect(colnames(motif2tf), tfs_use)
motif2tf <- motif2tf[, tfs_use]

rm(srt)
gc()

log_message('Fitting models.', verbose=verbose)
names(genes_use) <- genes_use
model_fits <- lapply_use(genes_use, function(g){

    gc()

    # Select peaks near gene
    gene_peaks <- as.logical(peaks2gene[g, ])
    if (sum(gene_peaks)==0){
        log_message('Warning: No peaks found near ', g, verbose=print_genes)
        return()
    }

    # Select peaks correlating with target gene expression
    peak_x <- peak_data[cluster_groups, gene_peaks, drop=F]
    peak_g_cor <- sparse_cov(peak_x, gene_data[, g, drop=F])$cor
    peak_g_cor[is.na(peak_g_cor)] <- 0
    peaks_use <- rownames(peak_g_cor)[abs(peak_g_cor[, 1]) > args$peak_cor]
    if (length(peaks_use)==0){
        log_message('Warning: No correlating peaks found for ', g, verbose=print_genes)
        return()
    }
    peak_x <- peak_x[, peaks_use, drop=F]
    peak_motifs <- peaks2motif[gene_peaks, , drop=F][peaks_use, , drop=F]

    # Select TFs with motifs in peaks
    gene_peak_tfs <- map(rownames(peak_motifs), function(p){
        x <- as.logical(peak_motifs[p, ])
        peak_tfs <- sparseMatrixStats::colMaxs(motif2tf[x, , drop=F])
        peak_tfs <- colnames(motif2tf)[as.logical(peak_tfs)]
        peak_tfs <- setdiff(peak_tfs, g)
        return(peak_tfs)
    })
    names(gene_peak_tfs) <- rownames(peak_motifs)

    # Check correlation of peaks with target gene
    gene_tfs <- purrr::reduce(gene_peak_tfs, union)
    tf_x <- gene_data[, gene_tfs, drop=F]
    tf_g_cor <- sparse_cov(tf_x, gene_data[, g, drop=F])$cor
    tf_g_cor[is.na(tf_g_cor)] <- 0
    tfs_use <- rownames(tf_g_cor)[abs(tf_g_cor[, 1]) > args$tf_cor]
    if (length(tfs_use)==0){
        log_message('Warning: No correlating TFs found for ', g, verbose=print_genes)
        return()
    }

    # Filter TFs and make formula string
    frml_string <- map(names(gene_peak_tfs), function(p){
        peak_tfs <- gene_peak_tfs[[p]]
        peak_tfs <- peak_tfs[peak_tfs%in%tfs_use]
        if (length(peak_tfs)==0){
            return()
        }
        peak_name <- str_replace_all(p, '-', '_')
        tf_name <- str_replace_all(peak_tfs, '-', '_')
        formula_str <- paste(
            paste(peak_name, args$inter_term, tf_name, sep=' '), collapse = ' + ')
        return(list(tfs=peak_tfs, frml=formula_str))
    })
    frml_string <- frml_string[!map_lgl(frml_string, is.null)]
    if (length(frml_string)==0){
        log_message('Warning: No valid peak:TF pairs found for ', g, verbose=print_genes)
        return()
    }

    target <- str_replace_all(g, '-', '_')
    model_frml <- as.formula(
        paste0(target, ' ~ ', paste0(map(frml_string, function(x) x$frml),  collapse=' + '))
    )

    # Get expression data
    nfeats <- sum(map_dbl(frml_string, function(x) length(x$tfs)))
    gene_tfs <- purrr::reduce(map(frml_string, function(x) x$tfs), union)
    gene_x <- gene_data[, union(g, gene_tfs), drop=F]
    model_mat <- as.data.frame(cbind(gene_x, peak_x))
    colnames(model_mat) <- str_replace_all(colnames(model_mat), '-', '_')

    log_message('Fitting model with ', nfeats, ' variables for ', g, verbose=print_genes)
    fit <- try(fit_model(
        model_frml, data=model_mat,
        family=gaussian, method=args$method, alpha=args$alpha
    ), silent=TRUE)
    if (any(class(fit)=='try-error')){
        log_message('Warning: Fitting model failed for ', g, verbose=print_genes)
        return()
    } else {
        fit$gof$nvariables <- nfeats
        return(fit)
    }
})
model_fits <- model_fits[!map_lgl(model_fits, is.null)]

coefs <- map_dfr(model_fits, function(x) x$coefs, .id='target')
gof <- map_dfr(model_fits, function(x) x$gof, .id='target')

log_message('Wrting output.', verbose=verbose)
coefs %>% write_tsv('coefs.tsv')
gof %>% write_tsv('gof.tsv')
