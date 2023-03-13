require(GenomeInfoDb)
require(EnsDb.Hsapiens.v86)
require(BSgenome.Hsapiens.UCSC.hg38)
require(JASPAR2020)
require(TFBSTools)
require(chromVAR)
require(GenomicRanges)
require(future.apply)
require(Signac)
library(fastmatch)


find_peaks_near_genes <- function(
    peaks, genes, distance = 200000, sep = c('-', '-'), only_tss = FALSE
){
    if (only_tss){
        genes <- resize(x = genes, width = 1, fix = 'start')
        genes_extended <- suppressWarnings(
            expr = Extend(
                genes, upstream = distance, downstream = distance
            )
        )
    } else {
        genes_extended <- suppressWarnings(
            expr = Extend(
                genes, upstream = distance
            )
        )
    }
    overlaps <- findOverlaps(
        query = peaks,
        subject = genes_extended,
        type = 'any',
        select = 'all'
    )
    hit_matrix <- sparseMatrix(
        i = queryHits(overlaps),
        j = subjectHits(overlaps),
        x = 1,
        dims = c(length(peaks), length(genes_extended))
    )
    rownames(hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
    colnames(hit_matrix) <- genes_extended$gene_name
    return(hit_matrix)
}


find_overlaps <- function(x, y, sep = c('-', '-')){
    overlaps <- findOverlaps(
        query = x,
        subject = y,
        type = 'any',
        select = 'all'
    )
    hit_matrix <- sparseMatrix(
        i = queryHits(overlaps),
        j = subjectHits(overlaps),
        x = 1,
        dims = c(length(x), length(y))
    )
    rownames(hit_matrix) <- GRangesToString(grange = x, sep = sep)
    colnames(hit_matrix) <- GRangesToString(grange = y, sep = sep)
    return(hit_matrix)
}

FindRegion <- function(
    object,
    region,
    sep = c("-", "-"),
    assay = NULL,
    extend.upstream = 0,
    extend.downstream = 0
){
    if (!is(object = region, class2 = "GRanges")) {
        # if separators are present in the string and we can convert the
        # start to a number, assume we're using genomic coordinates
        if (all(sapply(X = sep, FUN = grepl, x = region))) {
            region <- StringToGRanges(regions = region, sep = sep)
        } else {
            region <- LookupGeneCoords(object = object, assay = assay, gene = region)
            if (is.null(x = region)) {
                stop("Gene not found")
            }
        }
    }
    region <- suppressWarnings(
        expr = Extend(
            x = region,
            upstream = extend.upstream,
            downstream = extend.downstream
        )
    )
  return(region)
}

get_ranges_df <- function(
    ranges,
    region
){
    peak.intersect <- subsetByOverlaps(x = ranges, ranges = region)
    peak.df <- as.data.frame(x = peak.intersect)
    start.pos <- start(x = region)
    end.pos <- end(x = region)
    chromosome <- seqnames(x = region)
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    return(as_tibble(peak.df))
}

RangesPlot <- function(
    ranges,
    region,
    group.by = NULL,
    color = "dimgrey"
) {

    # subset to covered range
    peak.intersect <- subsetByOverlaps(x = ranges, ranges = region)
    peak.df <- as.data.frame(x = peak.intersect)
    start.pos <- start(x = region)
    end.pos <- end(x = region)
    chromosome <- seqnames(x = region)

    if (nrow(x = peak.df) > 0) {
        if (!is.null(x = group.by)) {
            if (!(group.by %in% colnames(x = peak.df))) {
                warning("Requested grouping variable not found")
                group.by <- NULL
        }
    }
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.plot <- ggplot(peak.df, aes_string(color=SetIfNull(x=group.by, y="color"))) +
        geom_segment(aes(x = start, y = 0, xend = end, yend = 0),
            size = 2,
            data = peak.df)
    } else {
        # no peaks present in region, make empty panel
        peak.plot <- ggplot(data = peak.df)
    }
    peak.plot <- peak.plot + theme_classic() +
        ylab(label = "Peaks") +
        theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank()) +
        xlab(label = paste0(chromosome, " position (bp)")) +
        xlim(c(start.pos, end.pos))

    if (is.null(x = group.by)) {
        # remove legend, change color
        peak.plot <- peak.plot +
            scale_color_manual(values = color) +
            theme(legend.position = "none")
    }
    return(peak.plot)
}

SetIfNull <- function(x, y) {
    if (is.null(x = x)) {
        return(y)
    } else {
        return(x)
    }
}


AnnotationPlot <- function(object, region) {
  annotation <- Annotation(object = object)
  if (is.null(x = annotation)) {
    return(NULL)
  }
  if (!inherits(x = region, what = 'GRanges')) {
    region <- StringToGRanges(regions = region)
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)

  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  genes.keep <- unique(x = annotation.subset$gene_name)
  annotation.subset <- annotation[
    fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
  ]

  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())

  } else {
    annotation.subset <- split(
      x = annotation.subset,
      f = annotation.subset$gene_name
    )
    p <- suppressWarnings(expr = suppressMessages(expr = ggbio::autoplot(
      object = annotation.subset,
      # GRangesFilter(value = region),
      fill = 'black',
      size = 1/2,
      color = 'black',
      names.expr = 'gene_name'
    )))
    p <- p@ggplot
    # extract y-axis limits and extend slightly so the label isn't covered
    y.limits <- ggplot_build(plot = p)$layout$panel_scales_y[[1]]$range$range
    p <- suppressMessages(p + ylim(y.limits[[1]], y.limits[[2]] + 0.5))
  }
  p <- p +
    theme_classic() +
    ylab('Genes') +
    xlab(label = paste0(chromosome, ' position (kb)')) +
    xlim(start.pos, end.pos) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  return(p)
}


DiffCoveragePlot <- function(
  object,
  region,
  group_by,
  ranges = NULL,
  annotation = TRUE,
  window = 1000,
  bins = 1000
){
    pc <- CoveragePlot(
        object,
        region = region,
        group.by = group_by,
        annotation = FALSE,
        peaks = FALSE,
        window = window,
    ) & no_legend() & theme(strip.text = element_blank())

    cov_df <- pc$data %>%
        mutate(cutpos=cut(position, bins, labels=F)) %>%
        group_by(cutpos, group) %>% summarize(cov=mean(coverage)) %>%
        pivot_wider(names_from='group', values_from=c('cov'), names_prefix='cov_') %>%
        mutate(coverage=cov_TRUE-cov_FALSE) %>%
        mutate(dir=factor(sign(coverage)))

    pd <- ggplot(cov_df, aes(cutpos, coverage)) +
        geom_hline(yintercept=0, color='darkgrey', size=0.1) +
        geom_area(stat='identity') +
        theme_void()

    pa <- AnnotationPlot(
        object,
        region = region
    ) + theme_void()


    if (!annotation){
        plt <- pd
    } else if (!is.null(ranges)){
        pr <- RangesPlot(
            ranges, region
        ) + theme_void() + no_legend()
        plt <- pd / pr / pa +
            plot_layout(heights=c(3,1,1)) & theme(plot.margin = margin(0,0,0,0))
    } else {
        plt <- pd / pa +
            plot_layout(heights=c(3,1)) & theme(plot.margin = margin(0,0,0,0))
    }
    return(plt)
}


#### Custom functions to get coverage from fragment files ####
get_peak_ranges <- function(object, assay='peaks', upstream=5000, downstream=5000){
    peak_ranges <- object[[assay]]@ranges
    peak_centers <- resize(peak_ranges, width=2, fix='center')
    peak_regions <- Extend(peak_centers, upstream = upstream, downstream = downstream)
    return(peak_regions)
}

get_coverage <- function(
  object, 
  regions, 
  group_by = NULL, 
  window = 500,
  bins = 500
){
  # Get coverage for all regions
  cov_plots <- map_par(1:length(regions), function(i){
    reg <- regions[i]
    p <- CoveragePlot(
      object, region=reg, annotation=F, peaks=F, group.by=group_by, window=window
    )
    return(p)
  })
  names(cov_plots) <- names(regions)
  cov_df <- map_dfr(cov_plots, ~.x$data, .id='region') %>%
    group_by(group, region) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, bins, labels=1:bins))) %>% 
    group_by(group, region, poscut) %>% 
    summarize(coverage=mean(coverage)) %>% 
    return()
}

# get_coverage <- function(object, ranges, binsize=10){
#     region_coverage_list <- Pando::map_par(1:length(ranges), function(i){
#         sparseMatrixStats::colMeans2(CutMatrix(object, region=ranges[i, ]))
#     }, parallel = T)
    
#     region_coverage <- Matrix(do.call(rbind, region_coverage_list), sparse=T)
    
#     if (binsize>1){
#         region_coverage <- t(Pando::aggregate_matrix(t(region_coverage), groups=cut(1:ncol(region_coverage), ncol(region_coverage)/10)))
#         colnames(region_coverage) <- 1:ncol(region_coverage)
#     }
#     return(region_coverage)
# }

# get_coverage_groups <- function(object, ranges, groups, binsize=10){
#     group_vec <- object[[groups]][,,drop=T]
#     region_coverage_list <- Pando::map_par(1:length(ranges), function(i){
#         Pando::aggregate_matrix(CutMatrix(object, region=ranges[i, ]), groups=group_vec)
#     }, parallel = T)
    
#     region_coverage <- Matrix(do.call(rbind, region_coverage_list), sparse=T)
    
#     if (binsize>1){
#         region_coverage <- t(Pando::aggregate_matrix(t(region_coverage), groups=cut(1:ncol(region_coverage), ncol(region_coverage)/10)))
#         colnames(region_coverage) <- 1:ncol(region_coverage)
#     }
    
#     region_coverage <- map(set_names(unique(group_vec)), function(x){
#         gcov <- region_coverage[rownames(region_coverage)==x, ]
#         rownames(gcov) <- NULL
#         return(gcov)
#     })
    
#     return(region_coverage)
# }



#### SIGNAC UTILS #### 
#### Not exported ####

#' @importFrom IRanges isDisjoint
NonOverlapping <- function(x, all.features) {
  # x is list of assays
  diff.features <- names(x = all.features[all.features < length(x = x)])
  if (length(x = diff.features) == 0) {
    return(TRUE)
  } else {
    diff.ranges <- StringToGRanges(regions = diff.features)
    return(isDisjoint(x = diff.ranges))
  }
}

#' @importFrom Matrix sparseMatrix
AddMissingCells <- function(x, cells) {
  # add columns with zeros for cells not in matrix
  missing.cells <- setdiff(x = cells, y = colnames(x = x))
  if (!(length(x = missing.cells) == 0)) {
    null.mat <- sparseMatrix(
      i = c(),
      j = c(),
      dims = c(nrow(x = x), length(x = missing.cells))
    )
    rownames(x = null.mat) <- rownames(x = x)
    colnames(x = null.mat) <- missing.cells
    x <- cbind(x, null.mat)
  }
  x <- x[, cells, drop = FALSE]
  return(x)
}

#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom Matrix Diagonal tcrossprod rowSums
AverageCountMatrix <- function(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL
) {
  assay = SetIfNull(x = assay, y = DefaultAssay(object = object))
  countmatrix <- GetAssayData(object = object[[assay]], slot = "counts")
  ident.matrix <- BinaryIdentMatrix(
    object = object,
    group.by = group.by,
    idents = idents
  )
  collapsed.counts <- tcrossprod(x = countmatrix, y = ident.matrix)
  avg.counts <- tcrossprod(
    x = collapsed.counts,
    y = Diagonal(x = 1 / rowSums(x = ident.matrix))
  )
  return(as.matrix(x = avg.counts))
}

# Create binary cell x class matrix of group membership
#' @importFrom Matrix sparseMatrix
BinaryIdentMatrix <- function(object, group.by = NULL, idents = NULL) {
  group.idents <- GetGroups(object = object, group.by = group.by, idents = idents)
  cell.idx <- seq_along(along.with = names(x = group.idents))
  unique.groups <- as.character(x = unique(x = group.idents))
  ident.idx <- seq_along(along.with = unique.groups)
  names(x = ident.idx) <- unique.groups
  ident.matrix <- sparseMatrix(
    i = ident.idx[as.character(x = group.idents)],
    j = cell.idx,
    x = 1
  )
  colnames(x = ident.matrix) <- names(x = group.idents)
  rownames(x = ident.matrix) <- unique.groups
  ident.matrix <- as(object = ident.matrix, Class = "dgCMatrix")
  return(ident.matrix)
}

# Calculate nCount and nFeature
#
# From Seurat
#
# @param object An Assay object
#
# @return A named list with nCount and nFeature
#
#' @importFrom Matrix colSums
#
CalcN <- function(object) {
  if (IsMatrixEmpty(x = GetAssayData(object = object, slot = "counts"))) {
    return(NULL)
  }
  return(list(
    nCount = colSums(x = object, slot = "counts"),
    nFeature = colSums(x = GetAssayData(object = object, slot = "counts") > 0)
  ))
}


#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

# Chunk GRanges
#
# Split a genomic ranges object into evenly sized chunks
# @param granges A GRanges object
# @param nchunk Number of chunks to split into
#
# @return Returns a list of GRanges objects
# @examples
# ChunkGRanges(blacklist_hg19, n = 10)
ChunkGRanges <- function(granges, nchunk) {
  if (length(x = granges) < nchunk) {
    nchunk <- length(x = granges)
  }
  chunksize <- as.integer(x = (length(granges) / nchunk))
  range.list <- sapply(X = seq_len(length.out = nchunk), FUN = function(x) {
    chunkupper <- (x * chunksize)
    if (x == 1) {
      chunklower <- 1
    } else {
      chunklower <- ((x - 1) * chunksize) + 1
    }
    if (x == nchunk) {
      chunkupper <- length(x = granges)
    }
    return(granges[chunklower:chunkupper])
  })
  return(range.list)
}

# Extract cell
#
# Extract cell barcode from list of tab delimited character
# vectors (output of \code{\link{scanTabix}})
#
# @param x List of character vectors
# @return Returns a string
#' @importFrom stringi stri_split_fixed
ExtractCell <- function(x) {
  if (length(x = x) == 0) {
    return(NULL)
  } else {
    x <- stri_split_fixed(str = x, pattern = "\t")
    n <- length(x = x)
    x <- unlist(x = x)
    return(unlist(x = x)[5 * (1:n) - 1])
  }
}

# Run groupCommand for the first n lines, convert the cell barcodes in the file
# to the cell names that appear in the fragment object, and subset the output to
# cells present in the fragment object
#
# Every cell in the fragment file will be present in the output dataframe. If
# the cell information is not set, every cell barcode that appears in the first
# n lines will be present.
#
# @param fragments A Fragment object
# @param n Number of lines to read from the beginning of the fragment file
# @param verbose Display messages
#
# @return Returns a data.frame
ExtractFragments <- function(fragments, n = NULL, verbose = TRUE) {
  fpath <- GetFragmentData(object = fragments, slot = "path")
  if (isRemote(x = fpath)) {
    stop("Remote fragment files not supported")
  }
  fpath <- normalizePath(path = fpath, mustWork = TRUE)
  cells <- GetFragmentData(object = fragments, slot = "cells")
  if (!is.null(x = cells)) {
    cells.use <- as.character(x = cells)
  } else {
    cells.use <- NULL
  }
  verbose <- as.logical(x = verbose)
  n <- SetIfNull(x = n, y = 0)
  n <- as.numeric(x = n)
  n <- round(x = n, digits = 0)
  counts <- groupCommand(
    fragments = fpath,
    some_whitelist_cells = cells.use,
    max_lines = n,
    verbose = verbose
  )
  # convert cell names
  if (!is.null(x = cells)) {
    # every cell will be present in the output, even if 0 counts
    converter <- names(x = cells)
    names(x = converter) <- cells
    counts$CB <- converter[counts$CB]
  }
  return(counts)
}

# convert region argument to genomic coordinates
# region can be a string, name of a gene, or GRanges object
FindRegion <- function(
  object,
  region,
  sep = c("-", "-"),
  assay = NULL,
  extend.upstream = 0,
  extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # first try to convert to coordinates, if not lookup gene
    region <- tryCatch(
      expr = suppressWarnings(
        expr = StringToGRanges(regions = region, sep = sep)
      ),
      error = function(x) {
        region <- LookupGeneCoords(
          object = object,
          assay = assay,
          gene = region
        )
        return(region)
      }
    )
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}

# GetReadsInRegion
#
# Extract reads for each cell within a given genomic region or set of regions
#
# @param cellmap A mapping of cell names in the fragment file to cell names in
# the Seurat object. Should be a named vector where each element is a cell name
# that appears in the fragment file and the name of each element is the
# name of the cell in the Seurat object.
# @param region A genomic region, specified as a string in the format
# 'chr:start-end'. Can be a vector of regions.
# @param tabix.file A TabixFile object.
# @param cells Cells to include. Default is all cells present in the object.
# @param verbose Display messages
# @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom Seurat Idents
#' @importFrom fastmatch fmatch
#
# @return Returns a data frame
GetReadsInRegion <- function(
  cellmap,
  region,
  tabix.file,
  cells = NULL,
  verbose = TRUE,
  ...
) {
  file.to.object <- names(x = cellmap)
  names(x = file.to.object) <- cellmap

  if (verbose) {
    message("Extracting reads in requested region")
  }
  if (!is(object = region, class2 = "GRanges")) {
    region <- StringToGRanges(regions = region, ...)
  }
  # remove regions that aren't in the fragment file
  common.seqlevels <- intersect(
    x = seqlevels(x = region),
    y = seqnamesTabix(file = tabix.file)
  )
  region <- keepSeqlevels(
    x = region,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  reads <- scanTabix(file = tabix.file, param = region)
  reads <- TabixOutputToDataFrame(reads = reads)
  reads <- reads[
    fmatch(x = reads$cell, table = cellmap, nomatch = 0L) > 0,
  ]
  # convert cell names to match names in object
  reads$cell <- file.to.object[reads$cell]
  if (!is.null(x = cells)) {
    reads <- reads[reads$cell %in% cells, ]
  }
  if (nrow(reads) == 0) {
    return(reads)
  }
  reads$length <- reads$end - reads$start
  return(reads)
}

# Get vector of cell names and associated identity
# @param object A Seurat object
# @param group.by Identity class to group cells by
# @param idents which identities to include
# @return Returns a named vector
#' @importFrom Seurat Idents
GetGroups <- function(
  object,
  group.by,
  idents
) {
  if (is.null(x = group.by)) {
    obj.groups <- Idents(object = object)
  } else {
    obj.md <- object[[group.by]]
    obj.groups <- obj.md[, 1]
    names(obj.groups) <- rownames(x = obj.md)
  }
  if (!is.null(idents)) {
    obj.groups <- obj.groups[obj.groups %in% idents]
  }
  return(obj.groups)
}

# Check if path is remote
# @param x path/s to check
isRemote <- function(x) {
  return(grepl(pattern = "^http|^ftp", x = x))
}

# row merge list of matrices
# @param mat.list list of sparse matrices
# @param new.rownames rownames to assign merged matrix
#' @importFrom Seurat RowMergeSparseMatrices
MergeMatrixParts <- function(mat.list, new.rownames) {
  # RowMergeSparseMatrices only exported in Seurat release Dec-2019 (3.1.2)
  merged.all <- mat.list[[1]]
  for (i in 2:length(x = mat.list)) {
    merged.all <- RowMergeSparseMatrices(
      mat1 = merged.all,
      mat2 = mat.list[[i]]
    )
  }
  # reorder rows to match genomic ranges
  merged.all <- merged.all[new.rownames, ]
  return(merged.all)
}

# Run GetReadsInRegion for a list of Fragment objects
# concatenate the output dataframes and return
# @param object A Seurat or ChromatinAssay object
# @param region Genomic region to extract fragments for
# @param fragment.list A list of Fragment objects. If NULL, pull them from the
# object
# @param assay Name of assay to use if supplying a Seurat object
#' @importFrom Seurat DefaultAssay
#' @importFrom Rsamtools TabixFile
#' @importFrom GenomeInfoDb keepSeqlevels
MultiGetReadsInRegion <- function(
  object,
  region,
  fragment.list = NULL,
  assay = NULL,
  ...
) {
  if (inherits(x = object, what = "Seurat")) {
    # pull the assay
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    object <- object[[assay]]
  }
  fragment.list <- SetIfNull(
    x = fragment.list,
    y = Fragments(object = object)
  )
  if (length(x = fragment.list) == 0) {
    # no fragments set
    stop("No fragment files found")
  }
  res <- data.frame()
  for (i in seq_along(along.with = fragment.list)) {
    tbx.path <- GetFragmentData(object = fragment.list[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragment.list[[i]], slot = "cells")
    tabix.file <- TabixFile(file = tbx.path)
    open(con = tabix.file)
    reads <- GetReadsInRegion(
      cellmap = cellmap,
      region = region,
      tabix.file = tabix.file,
      ...
    )
    res <- rbind(res, reads)
    close(con = tabix.file)
  }
  return(res)
}

# Generate matrix of integration sites
#
# Generates a cell-by-position matrix of Tn5 integration sites
# centered on a given region (usually a DNA sequence motif). This
# matrix can be used for downstream footprinting analysis.
#
# @param cellmap A mapping of cell names in the fragment file to cell names in
# the Seurat object. Should be a named vector where each element is a cell name
# that appears in the fragment file and the name of each element is the
# name of the cell in the Seurat object.
# @param region A set of GRanges containing the regions of interest
# @param cells Which cells to include in the matrix. If NULL, use all cells in
# the cellmap
# @param tabix.file A \code{\link[Rsamtools]{TabixFile}} object.
# @param verbose Display messages
#' @importFrom Matrix sparseMatrix
#' @importFrom Rsamtools TabixFile
#' @importMethodsFrom GenomicRanges width start end
# @return Returns a sparse matrix
SingleFileCutMatrix <- function(
  cellmap,
  region,
  cells = NULL,
  tabix.file,
  verbose = TRUE
) {
  # if multiple regions supplied, must be the same width
  cells <- SetIfNull(x = cells, y = names(x = cellmap))
  if (length(x = region) == 0) {
    return(NULL)
  }
  fragments <- GetReadsInRegion(
    region = region,
    cellmap = cellmap,
    cells = cells,
    tabix.file = tabix.file,
    verbose = verbose
  )
  start.lookup <- start(x = region)
  names(start.lookup) <- seq_along(region)
  # if there are no reads in the region
  # create an empty matrix of the correct dimension
  if (nrow(x = fragments) == 0) {
    cut.matrix <- sparseMatrix(
      i = NULL,
      j = NULL,
      dims = c(length(x = cells), width(x = region)[[1]])
    )
  } else {
    fragstarts <- start.lookup[fragments$ident] + 1
    cut.df <- data.frame(
      position = c(fragments$start, fragments$end) - fragstarts,
      cell = c(fragments$cell, fragments$cell),
      stringsAsFactors = FALSE
    )
    cut.df <- cut.df[
      (cut.df$position > 0) & (cut.df$position <= width(x = region)[[1]]),
      ]
    cell.vector <- seq_along(along.with = cells)
    names(x = cell.vector) <- cells
    cell.matrix.info <- cell.vector[cut.df$cell]
    cut.matrix <- sparseMatrix(
      i = cell.matrix.info,
      j = cut.df$position,
      x = 1,
      dims = c(length(x = cells), width(x = region)[[1]])
    )
  }
  rownames(x = cut.matrix) <- cells
  colnames(x = cut.matrix) <- seq_len(width(x = region)[[1]])
  return(cut.matrix)
}

# Generate matrix of integration sites
#
# Generates a cell-by-position matrix of Tn5 integration sites.
#
# @param object A Seurat object
# @param region A GRanges object containing the region of interest
# @param assay A name of assay to use. Must be a \code{\link{ChromatinAssay}}
# containing a list of \code{\link{Fragment}} objects.
# @param cells Which cells to include in the matrix. If NULL (default), use all
# cells in the object
# @param group.by Name of grouping variable to use
# @param verbose Display messages
# @return Returns a sparse matrix
#' @importFrom Seurat DefaultAssay
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom GenomeInfoDb keepSeqlevels
CutMatrix <- function(
  object,
  region,
  group.by = NULL,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  # run SingleFileCutMatrix for each fragment file and combine results
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  cells <- SetIfNull(x = cells, y = colnames(x = object))
  fragments <- Fragments(object = object[[assay]])
  res <- list()
  for (i in seq_along(along.with = fragments)) {
    fragment.path <- GetFragmentData(object = fragments[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragments[[i]], slot = "cells")
    tabix.file <- Rsamtools::TabixFile(file = fragment.path)
    open(con = tabix.file)
    # remove regions that aren't in the fragment file
    seqnames.in.both <- intersect(
      x = seqnames(x = region),
      y = Rsamtools::seqnamesTabix(file = tabix.file)
    )
    region <- GenomeInfoDb::keepSeqlevels(
      x = region,
      value = seqnames.in.both,
      pruning.mode = "coarse"
    )
    if (length(x = region) != 0) {
      cm <- SingleFileCutMatrix(
        region = region,
        cellmap = cellmap,
        tabix.file = tabix.file,
        cells = cells,
        verbose = FALSE
      )
      res[[i]] <- cm
    }
    close(con = tabix.file)
  }
  res <- Reduce(f = `+`, x = res)
  return(res)
}

# Generate cut matrix for many regions
#
# Run CutMatrix on multiple regions and add them together.
# Assumes regions are pre-aligned.
#
# @param object A Seurat object
# @param regions A set of GRanges
# @param group.by Name of grouping variable to use
# @param fragments A list of Fragment objects
# @param assay Name of the assay to use
# @param cells Vector of cells to include
# @param verbose Display messages
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom Seurat DefaultAssay
#' @importFrom GenomeInfoDb keepSeqlevels
MultiRegionCutMatrix <- function(
  object,
  regions,
  group.by = NULL,
  fragments = NULL,
  assay = NULL,
  cells = NULL,
  verbose = FALSE
) {
  if (inherits(x = object, what = "Seurat")) {
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    object <- object[[assay]]
  }
  fragments <- SetIfNull(x = fragments, y = Fragments(object = object))
  res <- list()
  if (length(x = fragments) == 0) {
    stop("No fragment files present in assay")
  }
  for (i in seq_along(along.with = fragments)) {
    frag.path <- GetFragmentData(object = fragments[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragments[[i]], slot = "cells")
    if (is.null(x = cellmap)) {
      cellmap <- colnames(x = object)
      names(x = cellmap) <- cellmap
    }
    tabix.file <- TabixFile(file = frag.path)
    open(con = tabix.file)
    # remove regions that aren't in the fragment file
    common.seqlevels <- intersect(
      x = seqlevels(x = regions),
      y = seqnamesTabix(file = tabix.file)
    )
    regions <- keepSeqlevels(
      x = regions,
      value = common.seqlevels,
      pruning.mode = "coarse"
    )
    cm <- SingleFileCutMatrix(
      cellmap = cellmap,
      tabix.file = tabix.file,
      region = regions,
      verbose = verbose
    )
    close(con = tabix.file)
    res[[i]] <- cm
  }
  # each matrix contains data for different cells at same positions
  # bind all matrices together
  res <- do.call(what = rbind, args = res)
  return(res)
}

# Create cut site pileup matrix
#
# For a set of aligned genomic ranges, find the total number of
# integration sites per cell per base.
#
# @param object A Seurat object
# @param regions A GRanges object
# @param assay Name of the assay to use
# @param cells Which cells to include. If NULL, use all cells
# @param verbose Display messages
#' @importMethodsFrom GenomicRanges strand
CreateRegionPileupMatrix <- function(
  object,
  regions,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  if (length(x = regions) == 0) {
    stop("No regions supplied")
  }
  # split into strands
  on_plus <- strand(x = regions) == "+" | strand(x = regions) == "*"
  plus.strand <- regions[on_plus, ]
  minus.strand <- regions[!on_plus, ]

  # get cut matrices for each strand
  if (verbose) {
    message("Finding + strand cut sites")
  }
  cut.matrix.plus <- MultiRegionCutMatrix(
    regions = plus.strand,
    object = object,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  if (verbose) {
    message("Finding - strand cut sites")
  }
  cut.matrix.minus <- MultiRegionCutMatrix(
    regions = minus.strand,
    object = object,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )

  # reverse minus strand and add together
  if (is.null(x = cut.matrix.plus)) {
    full.matrix <- cut.matrix.minus[, rev(x = colnames(x = cut.matrix.minus))]
  } else if (is.null(x = cut.matrix.minus)) {
    full.matrix <- cut.matrix.plus
  } else {
    full.matrix <- cut.matrix.plus + cut.matrix.minus[, rev(
      x = colnames(x = cut.matrix.minus)
    )]
  }
  # rename so 0 is center
  region.width <- width(x = regions)[[1]]
  midpoint <- round(x = (region.width / 2))
  colnames(full.matrix) <- seq_len(length.out = region.width) - midpoint
  return(full.matrix)
}

# Apply function to integration sites per base per group
#
# Perform colSums on a cut matrix with cells in the rows
# and position in the columns, for each group of cells
# separately.
#
# @param mat A cut matrix. See \code{\link{CutMatrix}}
# @param groups A vector of group identities, with the name
# of each element in the vector set to the cell name.
# @param fun Function to apply to each group of cells.
# For example, colSums or colMeans.
# @param group.scale.factors Scaling factor for each group. Should
# be computed using the number of cells in the group and the average number of
# counts in the group.
# @param normalize Perform sequencing depth and cell count normalization
# @param scale.factor Scaling factor to use. If NULL (default), will use the
# median normalization factor for all the groups.
ApplyMatrixByGroup <- function(
  mat,
  groups,
  fun,
  normalize = TRUE,
  group.scale.factors = NULL,
  scale.factor = NULL
) {
  if (normalize) {
    if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
      stop("If normalizing counts, supply group scale factors")
    }
  }
  all.groups <- as.character(x = unique(x = groups))
  if (any(is.na(x = groups))) {
    all.groups <- c(all.groups, NA)
  }
  ngroup <- length(x = all.groups)
  npos <- ncol(x = mat)

  group <- unlist(
    x = lapply(X = all.groups, FUN = function(x) rep(x, npos))
  )
  position <- rep(x = as.numeric(x = colnames(x = mat)), ngroup)
  count <- vector(mode = "numeric", length = npos * ngroup)

  for (i in seq_along(along.with = all.groups)) {
    grp <- all.groups[[i]]
    if (is.na(x = grp)) {
      pos.cells <- names(x = groups)[is.na(x = groups)]
    } else {
      pos.cells <- names(x = groups)[groups == all.groups[[i]]]
    }
    if (length(x = pos.cells) > 1) {
      totals <- fun(x = mat[pos.cells, ])
    } else {
      totals <- mat[pos.cells, ]
    }
    count[((i - 1) * npos + 1):((i * npos))] <- totals
  }

  # construct dataframe
  coverages <- data.frame(
    "group" = group, "position" = position, "count" = count,
    stringsAsFactors = FALSE
  )

  if (normalize) {
    scale.factor <- SetIfNull(
      x = scale.factor, y = median(x = group.scale.factors)
    )
    coverages$norm.value <- coverages$count /
      group.scale.factors[coverages$group] * scale.factor
  } else {
    coverages$norm.value <- coverages$count
  }
  return(coverages)
}

# TabixOutputToDataFrame
#
# Create a single dataframe from list of character vectors
#
# @param reads List of character vectors (the output of \code{\link{scanTabix}})
# @param record.ident Add a column recording which region the reads overlapped
# with
#' @importFrom stringi stri_split_fixed
#' @importFrom S4Vectors elementNROWS
TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  if (record.ident) {
    nrep <- elementNROWS(x = reads)
  }
  reads <- unlist(x = reads, use.names = FALSE)
  if (length(x = reads) == 0) {
    df <- data.frame(
      "chr" = "",
      "start" = "",
      "end" = "",
      "cell" = "",
      "count" = ""
    )
    df <- df[-1, ]
    return(df)
  }
  reads <- stringi::stri_split_fixed(str = reads, pattern = "\t")
  n <- length(x = reads[[1]])
  unlisted <- unlist(x = reads)
  e1 <- unlisted[n * (seq_along(along.with = reads)) - (n - 1)]
  e2 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 2)])
  e3 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 3)])
  e4 <- unlisted[n * (seq_along(along.with = reads)) - (n - 4)]
  e5 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 5)])
  df <- data.frame(
    "chr" = e1,
    "start" = e2,
    "end" = e3,
    "cell" = e4,
    "count" = e5,
    stringsAsFactors = FALSE,
    check.rows = FALSE,
    check.names = FALSE
  )
  if (record.ident) {
    df$ident <- rep(x = seq_along(along.with = nrep), nrep)
  }
  return(df)
}

# Extract delimiter information from a string.
#
# From Seurat
#
# Parses a string (usually a cell name) and extracts fields based on a delimiter
#
# @param string String to parse.
# @param field Integer(s) indicating which field(s) to extract. Can be a vector
# multiple numbers.
# @param delim Delimiter to use, set to underscore by default.
#
# @return A new string, that parses out the requested fields, and (if multiple),
# rejoins them with the same delimiter
#' @importFrom stringi stri_split_fixed
#
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(
    x = unlist(x = stri_split_fixed(
      str = as.character(x = field), pattern = ",")
    )
  )
  if (length(x = fields) == 1) {
    return(stri_split_fixed(str = string, pattern = delim)[[1]][field])
  }
  return(paste(
    stri_split_fixed(str = string, pattern = delim)[[1]][fields],
    collapse = delim))
}

# Set a default value if an object is null
#
# @param x An object to set if it's null
# @param y The value to provide if x is null
# @return Returns y if x is null, otherwise returns x.
SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

# Check if a matrix is empty
#
# From Seurat
#
# Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#
# @param x A matrix
#
# @return Whether or not \code{x} is empty
#
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

# Find matrix indices corresponding to overlapping genomic ranges
# @param assay.list A list of ChromatinAssay objects
# @param all.ranges Combined genomic ranges for all objects. This should be the
# set of ranges that \code{reduce} was run on to get \code{reduced.ranges}
# @param reduced.ranges A set of reduced genomic ranges containing the rev.map
# information
GetRowsToMerge <- function(assay.list, all.ranges, reduced.ranges) {
  revmap <- as.vector(x = reduced.ranges$revmap)

  # get indices of ranges that changed
  revmap.lengths <- sapply(X = revmap, FUN = length)
  changed.ranges <- which(x = revmap.lengths > 1)
  grange.string <- GRangesToString(grange = reduced.ranges[changed.ranges])

  # preallocate
  offsets <- list()
  results <- list()
  matrix.indices <- vector(
    mode = "numeric",
    length = length(x = changed.ranges) * 2
  )
  granges <- vector(
    mode = "character",
    length = length(x = changed.ranges) * 2
  )
  for (i in seq_along(along.with = assay.list)) {
    indices <- which(x = all.ranges$dataset == i)
    offsets[[i]] <- min(indices) - 1
    offsets[[i]][[2]] <- max(indices) + 1
    results[['matrix']][[i]] <- matrix.indices
    results[['grange']][[i]] <- granges
  }

  # find sets of ranges for each dataset
  counter <- vector(mode = "numeric", length = length(x = assay.list))
  for (x in seq_along(along.with = changed.ranges)) {
    idx <- changed.ranges[[x]]
    all.assay <- revmap[[idx]]
    for (i in seq_along(along.with = assay.list)) {
      this.assay <- all.assay[
        (all.assay > offsets[[i]][1]) & (all.assay < offsets[[i]][2])
      ]
      mat.idx <- this.assay - offsets[[i]][1]
      mat.idx <- mat.idx[mat.idx < offsets[[i]][2] & mat.idx > 0]
      for (y in seq_along(along.with = mat.idx)) {
        counter[i] <- counter[i] + 1
        results[['matrix']][[i]][[counter[i]]] <- mat.idx[[y]]
        results[['grange']][[i]][[counter[i]]] <- grange.string[[x]]
      }
    }
  }
  # remove trailing extra values in each vector
  for (i in seq_along(along.with = assay.list)) {
    results$matrix[[i]] <- results$matrix[[i]][1:counter[i]]
    results$grange[[i]] <- results$grange[[i]][1:counter[i]]
  }
  return(results)
}

# Merge rows of count matrices with overlapping genomic ranges
# @param mergeinfo The output of GetRowsToMerge: a list of matrix indices
#  and matrix rownames to be merged for each assay
# @param assay.list List of assays
# @param verbose Display messages
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix rowSums
#' @importMethodsFrom Matrix t
MergeOverlappingRows <- function(
  mergeinfo,
  assay.list,
  slot = "counts",
  verbose = TRUE
) {
  merge.counts <- list()
  for (i in seq_along(along.with = assay.list)) {
    # get count matrix
    counts <- GetAssayData(object = assay.list[[i]], slot = slot)

    if (nrow(x = counts) == 0) {
      # no counts, only data
      # skip row merge and return empty counts matrices
      merge.counts <- lapply(
        X = seq_along(along.with = assay.list),
        FUN = matrix,
        nrow = 0,
        ncol = 0
      )
      return(merge.counts)
    }

    # transpose for faster access since matrix is column major
    counts <- t(x = counts)

    # get rows to merge
    mrows <- mergeinfo$matrix[[i]]
    new.rownames <- mergeinfo$grange[[i]]
    nrep <- rle(x = new.rownames)

    # allocate
    todelete <- c()
    newmat <- vector(
      mode = "list",
      length = length(new.rownames)
    )
    newmat.names <- vector(
      mode = "character",
      length = length(x = new.rownames)
    )
    x <- 1  # row index for matrix
    y <- 1  # counter for list index
    if (verbose & length(x = nrep$lengths) > 1) {
      pb <- txtProgressBar(
        min = 1,
        max = length(x = nrep$lengths),
        style = 3,
        file = stderr()
      )
    }
    to.rename.idx <- vector(
      mode = "numeric", length = length(x = nrep$lengths)
    )
    to.rename.names <- vector(
      mode = "character", length = length(x = nrep$lengths)
    )
    idx.counter <- 0
    for (j in seq_along(along.with = nrep$lengths)) {
      rowrun <- nrep$lengths[[j]]
      new.feature.name <- nrep$values[[j]]
      index.range <- x:(x + rowrun - 1)
      matrix.index <- mrows[index.range]
      if (rowrun < 2) {
        idx.counter <- idx.counter + 1
        # no merge needed, just rename row in-place
        # store row indices and names to do the change in one step at the end
        to.rename.idx[idx.counter] <- matrix.index
        to.rename.names[idx.counter] <- new.feature.name
      } else {
        # merge multiple rows and add to list
        newmat[[y]] <- rowSums(x = counts[, matrix.index])
        # mark merged row for deletion
        todelete <- c(todelete, matrix.index)
        # add row names
        newmat.names[y] <- new.feature.name
        y <- y + 1
      }
      if (verbose & length(x = nrep$lengths) > 1) {
        setTxtProgressBar(pb = pb, value = j)
      }
      x <- x + rowrun
    }
    # remove extra elements in vectors
    to.rename.idx <- to.rename.idx[1:idx.counter]
    to.rename.names <- to.rename.names[1:idx.counter]
    newmat <- newmat[1:(y - 1)]
    newmat.names <- newmat.names[1:(y - 1)]

    # transpose back
    counts <- t(x = counts)

    # rename matrix rows that weren't merged
    rownames(counts)[to.rename.idx] <- to.rename.names

    if (y == 1) {
      # no rows were merged, can return counts
      merge.counts[[i]] <- counts
    } else if (y == 2) {
      # only one element
      tomerge <- matrix(data = newmat[[1]], nrow = 1)
      colnames(x = tomerge) <- names(x = newmat[[1]])
      rownames(x = tomerge) <- newmat.names
      tomerge <- tomerge[, colnames(x = counts)]
      counts <- rbind(counts, tomerge)
      merge.counts[[i]] <- counts
    } else {
      # construct sparse matrix
      if (verbose) {
        message("\nBinding matrix rows")
      }
      merged.mat <- Reduce(f = rbind, x = newmat)
      rownames(merged.mat) <- newmat.names
      merged.mat <- as(object = merged.mat, Class = "dgCMatrix")

      # remove rows from count matrix that were merged
      mat.rows <- seq_len(length.out = nrow(x = counts))
      tokeep <- setdiff(mat.rows, todelete)
      counts <- counts[tokeep, ]

      # add new merged rows to counts
      counts <- rbind(counts, merged.mat)
      merge.counts[[i]] <- counts
    }
  }
  return(merge.counts)
}

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors elementNROWS
PartialMatrix <- function(tabix, regions, sep = c("-", "-"), cells = NULL) {
  # construct sparse matrix for one set of regions
  # names of the cells vector can be ignored here, conversion is handled in
  # the parent functions
  open(con = tabix)
  cells.in.regions <- GetCellsInRegion(
    tabix = tabix,
    region = regions,
    cells = cells
  )
  close(con = tabix)
  gc(verbose = FALSE)
  nrep <- elementNROWS(x = cells.in.regions)
  if (all(nrep == 0) & !is.null(x = cells)) {
    # no fragments
    # zero for all requested cells
    featmat <- sparseMatrix(
      dims = c(length(x = regions), length(x = cells)),
      i = NULL,
      j = NULL
    )
    rownames(x = featmat) <- GRangesToString(grange = regions)
    colnames(x = featmat) <- cells
    featmat <- as(object = featmat, Class = "dgCMatrix")
    return(featmat)
  } else if (all(nrep == 0)) {
    # no fragments, no cells requested
    # create empty matrix
    featmat <- sparseMatrix(
      dims = c(length(x = regions), 0),
      i = NULL,
      j = NULL
    )
    rownames(x = featmat) <- GRangesToString(grange = regions)
    featmat <- as(object = featmat, Class = "dgCMatrix")
    return(featmat)
  } else {
    # fragments detected
    if (is.null(x = cells)) {
      all.cells <- unique(x = unlist(x = cells.in.regions))
      cell.lookup <- seq_along(along.with = all.cells)
      names(x = cell.lookup) <- all.cells
    } else {
      cell.lookup <- seq_along(along.with = cells)
      names(cell.lookup) <- cells
    }
    # convert cell name to integer
    cells.in.regions <- unlist(x = cells.in.regions)
    cells.in.regions <- unname(obj = cell.lookup[cells.in.regions])
    all.features <- GRangesToString(grange = regions, sep = sep)
    feature.vec <- rep(x = seq_along(along.with = all.features), nrep)
    featmat <- sparseMatrix(
      i = feature.vec,
      j = cells.in.regions,
      x = rep(x = 1, length(x = cells.in.regions))
    )
    featmat <- as(Class = "dgCMatrix", object = featmat)
    rownames(x = featmat) <- all.features[1:max(feature.vec)]
    colnames(x = featmat) <- names(x = cell.lookup)[1:max(cells.in.regions)]
    # add zero columns for missing cells
    if (!is.null(x = cells)) {
      featmat <- AddMissingCells(x = featmat, cells = cells)
    }
    # add zero rows for missing features
    missing.features <- all.features[!(all.features %in% rownames(x = featmat))]
    if (length(x = missing.features) > 0) {
      null.mat <- sparseMatrix(
        i = c(),
        j = c(),
        dims = c(length(x = missing.features), ncol(x = featmat))
      )
      rownames(x = null.mat) <- missing.features
      null.mat <- as(object = null.mat, Class = "dgCMatrix")
      featmat <- rbind(featmat, null.mat)
    }
    return(featmat)
  }
}

# Convert PFMMatrix to
# @param x A PFMatrix
PFMatrixToList <- function(x) {
  if (!requireNamespace("TFBSTools", quietly = TRUE)) {
    stop("Please install TFBSTools.
         https://www.bioconductor.org/packages/TFBSTools/")
  }
  position.matrix <- TFBSTools::Matrix(x = x)
  name.use <- TFBSTools::name(x = x)
  return(list("matrix" = position.matrix, "name" = name.use))
}

#' @importFrom Matrix rowMeans rowSums
SparseRowVar <- function(x) {
  return(rowSums(x = (x - rowMeans(x = x)) ^ 2) / (dim(x = x)[2] - 1))
}

#' @importMethodsFrom Matrix t
SparseColVar <- function(x) {
  return(SparseRowVar(x = t(x = x)))
}


CollapseToLongestTranscript <- function(ranges) {
    range.df <- as.data.table(x = ranges)
    range.df$strand <- as.character(x = range.df$strand)
    range.df$strand <- ifelse(
        test = range.df$strand == "*",
        yes = "+",
        no = range.df$strand
    )
    collapsed <- range.df[
        , .(unique(seqnames),
            min(start),
            max(end),
            strand[[1]],
            gene_biotype[[1]],
            gene_name[[1]]),
        "gene_id"
    ]
    colnames(x = collapsed) <- c(
        "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
    )
    collapsed$gene_name <- make.unique(names = collapsed$gene_name)
    gene.ranges <- makeGRangesFromDataFrame(
        df = collapsed,
        keep.extra.columns = TRUE
    )
    return(gene.ranges)
}
