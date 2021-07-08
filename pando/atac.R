require(ChIPseeker)
require(GenomeInfoDb)
require(EnsDb.Hsapiens.v86)
require(BSgenome.Hsapiens.UCSC.hg38)
require(JASPAR2020)
require(TFBSTools)
require(chromVAR)
require(GenomicRanges)
require(future.apply)
require(Signac)
require(fastmatch)


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
      GRangesFilter(value = region),
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
