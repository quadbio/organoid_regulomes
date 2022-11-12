source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')

setwd('~/projects/early/')


#### FUNC ####
PlotGRNFeatures <- function(
    object, region, cons_regions, cre_regions, exon_ranges, assay='peaks', extend.upstream=50000, extend.downstream=0, group.by=NULL, cov_colors=many
){
    
    region <- FindRegion(
        object = object,
        region = region,
        assay = assay,
        extend.upstream = extend.upstream,
        extend.downstream = extend.downstream
    )
    
    pp <- PeakPlot(
        object = object,
        region = region
    ) + theme_void() + no_legend() +
        ggtitle(label='Peaks')
    
    pcre <- RangesPlot(
        ranges = cre_regions,
        region = region,
        color = '#004B25'
    ) + theme_void() + no_legend() +
        ggtitle(label='cis Regulatory Elements (cRE, SCREEN)')
    
    pcons <- RangesPlot(
        ranges = cons_regions,
        region = region,
        color = '#253894'
    ) + theme_void() + no_legend() +
        ggtitle(label='Conserved Elements (phastCons 20Mammals)')
    
    pex <- RangesPlot(
        ranges = exon_ranges,
        region = region,
        color = '#005667'
    ) + theme_void() + no_legend() +
        ggtitle(label='Exons')
    
    pgene <- AnnotationPlot(object, region=region) +
        theme_void()
    
    pout <- pp / pcre / pcons / pex / pgene + plot_layout(heights=c(1,1,1,1,2))
    return(pout)
}

#### Read stuff ####
rnatac <- read_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2_srt.rds')
# rnatac %>% write_rds('data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2_srt.rds')

atac <- read_rds('data/ATAC/ATAC_all_merged/ATAC_all_merged_v3.1_srt.rds')
rna <- read_rds('data/RNA/RNA_all_merged/RNA_all_merged_v2_srt.rds')
# motif_df$tf %>% write('~/resources/JASPAR2020_tfs.txt')
tfs <- read_tsv('~/resources/DB/animal_tfdb/tf_human.tsv')
jaspar_motifs <- read_rds('~/resources/DB/JASPAR/JASPAR2020_hs_all_collections_motifs.rds')

jaspar_df <- read_tsv('~/resources/mappings/motif2tf_JASPAR2020_all.tsv')
cbp_pfml <- read_rds('~/resources/DB/CIS-BP/CISPB_missing_motifs.rds')
cbp_df <- read_tsv('~/resources/DB/CIS-BP/motif2tf_CISPB_missing.rds')
needle_df <- read_tsv('data/tf_matches/motif2tf_early_var4k_missing.tsv')
all_df <- read_tsv('data/tf_matches/motif2tf_early_all.tsv')



#### Plots on motif selection ####
plot_df <- all_df %>% 
    distinct(origin, tf, exp_evidence) %>% 
    mutate(origin=factor(origin, levels=c('JASPAR2020', 'CIS-BP', 'SEQ_SIMILARITY')))
ggplot(plot_df, aes(origin, fill=exp_evidence)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('#424949','darkgrey', 'darkgrey')) +
    theme_rangeframe() + scale_axis_rangeframe() + 
    article_text() + no_legend() +
    scale_y_continuous(breaks=c(0, 500, 1000), limits=c(0,1000)) +
    rotate_x_text(40) +
    theme(
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    labs(y='Count')
ggsave('plots/paper/supp_grn/motif_tf_origin_bar.pdf', width=3, height=4, unit='cm')


ggplot(needle_df, aes(prc_similarity)) +
    geom_histogram(color='black', fill='darkgrey', size=0.2) +
    theme_rangeframe() + scale_axis_rangeframe() + 
    article_text() +
    labs(y='Count', x='Sequence similarity')
ggsave('plots/paper/supp_grn/motif_tf_seq_sim_hist.pdf', width=5, height=4, unit='cm')  
           


## Plot motifs for several families

# CP2-related factors
motifs_plot <- needle_df %>% filter(tf=='GRHL3') %>% pull(motif)
p1 <- MotifPlot(rnatac, motifs_plot, assay='peaks') &
    article_text() 

# bHLH
motifs_plot <- needle_df %>% filter(tf=='HES3') %>% pull(motif)
p2 <- MotifPlot(rnatac, motifs_plot, assay='peaks') &
    article_text() 

# Homeobox
motifs_plot <- needle_df %>% filter(tf=='HOPX') %>% pull(motif)
p3 <- MotifPlot(rnatac, motifs_plot, assay='peaks') &
    article_text() 

# Paired-related
motifs_plot <- needle_df %>% filter(tf=='GATAD2A') %>% pull(motif)
p4 <- MotifPlot(rnatac, motifs_plot, assay='peaks') &
    article_text() 

p1 / p2 / p3 & 
    theme_void()

ggsave('plots/paper/supp_grn/motif_tf_seq_sim_motif_plots.pdf', width=6, height=3.8, unit='cm')  



#### Plots on region filtering ####
load('data/phastConsElements20Mammals.UCSC.hg38.RData')
load('data/SCREEN.ccRE.UCSC.hg38.RData')
peak_subsets <- read_rds('data/RNA_ATAC/RNA_ATAC_peak_subsets_v2.rds')
peak_ranges <- rnatac@assays$peaks@ranges


exon_ranges <- gene_annot[gene_annot$type=='exon', ]
names(exon_ranges@ranges) <- NULL
exon_ranges <- intersect(exon_ranges, exon_ranges)
exon_ranges <- GRanges(
    seqnames = exon_ranges@seqnames,
    ranges = exon_ranges@ranges
)


# Compare totl peaks / peaks in exons / peaks in cRE / peaks in cons
peak_cons_overlaps <- findOverlaps(phastConsElements20Mammals.UCSC.hg38, peak_ranges)
peak_cons_matches <- subjectHits(peak_cons_overlaps)
peak_cons_matches %>% unique %>% length

peak_ex_overlaps <- findOverlaps(exon_ranges, peak_ranges)
peak_ex_matches <- subjectHits(peak_ex_overlaps)
peak_ex_matches %>% unique %>% length

peak_cre_overlaps <- findOverlaps(SCREEN.ccRE.UCSC.hg38, peak_ranges)
peak_cre_matches <- subjectHits(peak_cre_overlaps)
peak_cre_matches %>% unique %>% length

peak_cons_cre_overlaps <- findOverlaps(
    IRanges::union(phastConsElements20Mammals.UCSC.hg38, SCREEN.ccRE.UCSC.hg38),
    peak_ranges
)
peak_cons_cre_matches <- subjectHits(peak_cons_cre_overlaps)
peak_cons_cre_matches %>% unique %>% length


peak_cons_and_cre_overlaps <- findOverlaps(
    IRanges::intersect(phastConsElements20Mammals.UCSC.hg38, SCREEN.ccRE.UCSC.hg38),
    peak_ranges
)
peak_cons_and_cre_matches <- subjectHits(peak_cons_and_cre_overlaps)
peak_cons_and_cre_matches %>% unique %>% length


plot_df <- tibble(
    peaks = length(peak_ranges),
    peak_cons = peak_cons_matches %>% unique %>% length,
    peak_ex = peak_ex_matches %>% unique %>% length,
    peak_cre = peak_cre_matches %>% unique %>% length,
    peak_cons_cre = peak_cons_cre_matches %>% unique %>% length,
    peak_cons_and_cre = peak_cons_and_cre_matches %>% unique %>% length
) %>% pivot_longer(everything()) %>% 
    arrange(value) %>% mutate(name=factor(name, levels=unique(.$name)))

ggplot(plot_df, aes(name, value/1000)) +
    geom_bar(stat='identity', color='black', fill='darkgrey', size=0.2) +
    coord_flip() 

ggplot(plot_df, aes(name, value/1000)) +
    geom_bar(stat='identity', color='black', fill='darkgrey', size=0.2) +
    geom_text(mapping=aes(label=value), hjust=1, nudge_y=-10, size=5/ggplot2::.pt) +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_y_continuous(limits=c(0, 600)) +
    theme(
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
    ) +
    coord_flip() +
    article_text() +
    no_margin() +
    labs(y='Count')
ggsave('plots/paper/supp_grn/regions_intersect_bar.pdf', width=4, height=3.5, unit='cm')



# Some tracks on region overlaps
goi <- 'DLX6'
extend_upstream <-  10000
extend_downstream <- 1000

atac$all <- 'all'
pcov <- CoveragePlot(
    object = atac,
    group.by = 'all',
    region = goi,
    extend.upstream = extend_upstream,
    extend.downstream = extend_downstream,
    peaks = FALSE,
    annotation = FALSE
) & scale_fill_manual(values='black') & theme_void() & no_legend() & theme(strip.text = element_blank())

preg <- PlotGRNFeatures(
    atac, cons_regions = phastConsElements20Mammals.UCSC.hg38, cre_regions = SCREEN.ccRE.UCSC.hg38, exon_ranges = exon_ranges,
    region = goi, extend.upstream = extend_upstream, extend.downstream = extend_downstream,
) & theme_void() & theme(title=element_text(size=5)) & no_legend()

pcov / preg + plot_layout(heights=c(1,2))
ggsave('plots/paper/supp_grn/regions_dlx6_tracks.pdf', width=4.6, height=3.6, units='cm')


goi <- 'EMX1'
extend_upstream <-  20000
extend_downstream <- 1000

atac$all <- 'all'
pcov <- CoveragePlot(
    object = atac,
    group.by = 'all',
    region = goi,
    extend.upstream = extend_upstream,
    extend.downstream = extend_downstream,
    peaks = FALSE,
    annotation = FALSE
) & scale_fill_manual(values='black') & theme_void() & no_legend() & theme(strip.text = element_blank())

preg <- PlotGRNFeatures(
    atac, cons_regions = phastConsElements20Mammals.UCSC.hg38, cre_regions = SCREEN.ccRE.UCSC.hg38, exon_ranges = exon_ranges,
    region = goi, extend.upstream = extend_upstream, extend.downstream = extend_downstream,
) & theme_void() & theme(title=element_text(size=5)) & no_legend()

pcov / preg + plot_layout(heights=c(1,2))
ggsave('plots/paper/supp_grn/regions_emx1_tracks.pdf', width=4.6, height=3.6, units='cm')




to_tf <- 'DLX5'

annot <- gene_annot[gene_annot$gene_name==to_tf]
strand(annot) <- '-'
annot_ex <- Extend(annot, upstream = 100000)
region <- GRanges(seqnames = as.character(seqnames(annot_ex))[[1]], ranges = IRanges(start = min(start(annot_ex)), end = max(end(annot_ex))))
# region <- resize(region, width = 30000, fix='start')


binding_regions <- modules %>% filter(target==to_tf) %>% pull(peaks) %>% str_split(';') %>% unlist() %>% StringToGRanges(sep=c('_', '_'))

binding_regions <- modules %>% filter(target==to_tf) %>% 
    separate_rows(peaks, sep=';') %>% 
    separate(peaks, sep='_', into=c('seqnames', 'start', 'end'), convert=T) %>% 
    filter(end>start(region), start<end(region)) %>% 
    group_by(start, end) %>% 
    mutate(padj_peak=min(padj)) %>% 
    mutate(padj_peak=pmax(padj, 1e-10)) %>% 
    mutate(pos_ratio=sum(sign(estimate)==1)/n()) %>% 
    ungroup()





#### Plots with fit metrics ####
coefs <- read_tsv('data/RNA_ATAC/grn/v5/vall_d100k_cons_cre_noex_glm/coefs.tsv')
gof <- read_tsv('data/RNA_ATAC/grn/v5/vall_d100k_cons_cre_noex_glm/gof.tsv')
modules <- read_tsv('data/RNA_ATAC/grn/networks/v5_d100k_cons_cre_noex_modules.tsv')

#### Reformat and filter ####

gof <- gof %>% dplyr::filter(dsq<=1, dsq>=0) %>% 
    mutate(nice=dsq>0.10&nvariables>10) 

models_use <- gof %>% filter(nice) %>% pull(target) %>% unique()
coefs_use <- coefs %>% 
    mutate(padj = p.adjust(pval, method = 'fdr')) %>% 
    filter(target%in%models_use) %>% 
    dplyr::filter(term!='(Intercept)') %>% 
    mutate(
        tf_ = str_replace(term, '(.+):(.+)', '\\1'),
        peak_ = str_replace(term, '(.+):(.+)', '\\2')
    ) %>% 
    mutate(
        tf = ifelse(str_detect(tf_, 'chr'), peak_, tf_),
        peak = ifelse(!str_detect(tf_, 'chr'), peak_, tf_)
    ) %>% dplyr::select(-peak_, -tf_) 



## GOF metrics
p1 <- ggplot(gof, aes(dsq, nvariables, alpha=nice)) +
    geom_pointdensity(size=0.01, shape=16) +
    geom_hline(yintercept=10, size=0.2, color='darkgrey', linetype='dashed') +
    geom_vline(xintercept=0.1, size=0.2, color='darkgrey', linetype='dashed') +
    scale_color_gradientn(colors=rev(reds(1.3))) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 1, 10, 100, 1000, 10000)) +
    scale_x_continuous(breaks=seq(0,1,0.2)) +
    article_text() +
    no_legend() +
    labs(x=expression('Explained variance'~(R**2)), y='# variables in model') +
    theme(
        plot.margin = unit(c(0,0,0,0), 'line'),
        strip.text = element_blank()
    )

p2 <- ggplot(gof, aes(dsq)) +
    geom_histogram(fill='darkgray', bins=20, color='black', size=0.2) +
    theme_void() +
    no_legend()

p3 <- ggplot(gof, aes(nvariables)) +
    geom_histogram(fill='darkgray', bins=20, color='black', size=0.2) +
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 1, 10, 100, 1000, 10000)) +
    theme_void() +
    coord_flip() +
    no_legend() 


layout <- '
AAAA#
BBBBC
BBBBC
'
p2 + p1 + p3 + plot_layout(design = layout) & no_margin()
ggsave('plots/paper/supp_grn/model_nvar_vs_r2_scatter.pdf', width=4.8, height=3.3, units='cm')



## Module sizes
plot_df <- modules %>% 
    distinct(target, npeaks)

p1 <- ggplot(plot_df, aes(1, npeaks)) +
    geom_violin(size=0.2, fill='darkgrey', color='black') +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
    labs(y='# peaks')


plot_df <- modules %>% 
    distinct(target, ntfs)

p2 <- ggplot(plot_df, aes(1, ntfs)) +
    geom_violin(size=0.2, fill='darkgrey', color='black') +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
    labs(y='# TFs')


plot_df <- modules %>% 
    distinct(tf, ngenes)

p3 <- ggplot(plot_df, aes(1, ngenes)) +
    geom_violin(size=0.2, fill='darkgrey', color='black') +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 1, 10, 100, 1000, 10000)) +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
    labs(y=expression('# genes'))

p1 | p2 | p3 & no_margin()
ggsave('plots/paper/supp_grn/model_module_specs_violin.pdf', width=5, height=3.3, units='cm')




#### Module activity feature plots ####

mod_plot <- c('POU5F1', 'PAX6', 'EMX1', 'NFIA', 'DLX5', 'SP8', 'LHX1', 'LMX1A')

rnatac@active.assay <- 'RNA'
p0 <- feature_plot(rnatac, features=mod_plot, order=T, ncol=length(mod_plot), pt.size=0.1) &
    scale_color_gradientn(colors=gyylgnbu()) & no_margin() & theme(text=element_blank())
rnatac@active.assay <- 'module_score'
p1 <- feature_plot(rnatac, features=mod_plot, order=T, ncol=length(mod_plot), pt.size=0.1) &
    scale_color_gradientn(colors=grad(pals::ocean.tempo, 0.5), na.value=pals::ocean.tempo(100)[0]) & no_margin() & theme(text=element_blank())
rnatac@active.assay <- 'module_neg_score'
p2 <- feature_plot(rnatac, features=mod_plot, order=T, ncol=length(mod_plot), pt.size=0.1) &
    scale_color_gradientn(colors=grad(pals::ocean.tempo, 0.5), na.value=pals::ocean.tempo(100)[0]) & no_margin() & theme(text=element_blank())
rnatac@active.assay <- 'module_peaks_chromvar'
p3 <- feature_plot(rnatac, features=mod_plot, order=T, ncol=length(mod_plot), pt.size=0.1) &
    scale_color_gradientn(colors=gyylorbr(1), na.value=gyylorbr(1)[0]) & no_margin() & theme(text=element_blank())
rnatac@active.assay <- 'module_neg_peaks_chromvar'
p4 <- feature_plot(rnatac, features=mod_plot, order=T, ncol=length(mod_plot), pt.size=0.1) &
    scale_color_gradientn(colors=gyylorbr(1), na.value=gyylorbr(1)[0]) & no_margin() & theme(text=element_blank())


p0 / p1 / p2 / p3 / p4
ggsave('plots/paper/supp_grn/modules_score_umap.png', width=12*4, height=6*4, unit='cm')



rnatac@active.assay <- 'RNA'
feature_plot(rnatac, features='DLX6', order=T, pt.size=0.1) &
    scale_color_gradientn(colors=gyylgnbu()) 

#### Trans probs on UMAP ####
umap_meta <- Reductions(rnatac, slot='umap')@cell.embeddings %>% 
  as_tibble(rownames='cell')

meta <- rnatac@meta.data %>% 
  as_tibble(rownames='cell') %>% 
  inner_join(umap_meta)

ggplot(meta, aes(UMAP_1, UMAP_2, color=to_ctx_ranks)) +
  geom_point(size=0.2) +
  theme_void() + no_legend() +
  scale_color_gradientn(colors=rev(pals::ocean.deep(100)))
 ggsave('plots/paper/supp_grn/trans_probs_ctx.png', width=5, height=4)
 
ggplot(meta, aes(UMAP_1, UMAP_2, color=to_ge_ranks)) +
  geom_point(size=0.2) +
  theme_void() + no_legend() +
  scale_color_gradientn(colors=rev(pals::ocean.deep(100)))
ggsave('plots/paper/supp_grn/trans_probs_ge.png', width=5, height=4)

ggplot(meta, aes(UMAP_1, UMAP_2, color=to_nt_ranks)) +
  geom_point(size=0.2) +
  theme_void() + no_legend() +
  scale_color_gradientn(colors=rev(pals::ocean.deep(100)))
ggsave('plots/paper/supp_grn/trans_probs_nt.png', width=5, height=4)





#### Predicted regulatory regions and tracks ####
atac_cov <- subset(atac, lineage%in%c('nt', 'telencephalon', 'early', 'ge', 'ctx'))
atac_cov$lineage <- factor(atac_cov$lineage, levels=c('early', 'nt', 'telencephalon', 'ge', 'ctx'))


#### POU5F1 ####

to_tf <- 'POU5F1'

annot <- gene_annot[gene_annot$gene_name==to_tf]
strand <- mode(as.character(strand(annot)))
region <- GRanges(
    seqnames = as.character(seqnames(annot))[[1]], 
    ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
    strand = strand
)
region <- resize(region, width=0, fix='start')
region <- Extend(region, upstream = 15000, downstream=5000)


binding_regions <- modules %>% filter(target==to_tf) %>% pull(peaks) %>% str_split(';') %>% unlist() %>% StringToGRanges(sep=c('_', '_'))

binding_regions <- modules %>% filter(target==to_tf) %>% 
    separate_rows(peaks, sep=';') %>% 
    separate(peaks, sep='_', into=c('seqnames', 'start', 'end'), convert=T) %>% 
    filter(end>start(region), start<end(region)) %>% 
    group_by(start, end) %>% 
    mutate(padj_peak=min(padj)) %>% 
    mutate(padj_peak=pmax(padj, 1e-10)) %>% 
    mutate(pos_ratio=sum(sign(estimate)==1)/n()) %>% 
    ungroup()


ggplot(binding_regions) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=4, max.overlaps = 100, nudge_y=-100, angle=90, force=5
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


binding_regions_label <- binding_regions %>% 
    filter(padj_peak<1e-4) %>% 
    group_by(start, end) %>% 
    top_n(3, estimate)

ggplot(binding_regions_label) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


plot_df <- binding_regions %>% 
    group_by(start, end) %>% 
    summarize(padj_peak=min(padj_peak))
p_mot <- ggplot(plot_df) +
    geom_text_repel(
        data = binding_regions_label,
        mapping = aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    ggnewscale::new_scale_color() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    xlim(c(start(region), end(region))) + no_legend()

p_cov <- CoveragePlot(
    atac_cov, 
    region = region, 
    group.by = 'lineage', 
    annotation = FALSE,
    peaks = FALSE
) & scale_fill_manual(values=lineage_colors) &
    theme_void() & no_legend() & theme(strip.text = element_blank())

p_an <- AnnotationPlot(
    atac_cov,
    region = region
) + theme_void() 

p_pe <- PeakPlot(
    atac_cov,
    region = region
) + theme_void() + no_legend()

p_cov / p_an / p_pe / p_mot + plot_layout(heights=c(5,1,1,3)) & theme(text = element_text(size=5))
ggsave('plots/paper/supp_grn/grn_reg_regions_pou5f1_tracks.pdf', width=4.4, height=3.7, unit='cm')



#### NEUROD6 ####

to_tf <- 'NEUROD6'
annot <- gene_annot[gene_annot$gene_name==to_tf]
strand <- mode(as.character(strand(annot)))
# strand <- '+'
region <- GRanges(
    seqnames = as.character(seqnames(annot))[[1]], 
    ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
    strand = strand
)
region <- resize(region, width=0, fix='start')
region <- Extend(region, upstream = 25000, downstream=1000)


binding_regions <- modules %>% filter(target==to_tf) %>% pull(peaks) %>% str_split(';') %>% unlist() %>% StringToGRanges(sep=c('_', '_'))

binding_regions <- modules %>% filter(target==to_tf) %>% 
    separate_rows(peaks, sep=';') %>% 
    separate(peaks, sep='_', into=c('seqnames', 'start', 'end'), convert=T) %>% 
    filter(end>start(region), start<end(region)) %>% 
    group_by(start, end) %>% 
    mutate(padj_peak=min(padj)) %>% 
    mutate(padj_peak=pmax(padj, 1e-10)) %>% 
    mutate(pos_ratio=sum(sign(estimate)==1)/n()) %>% 
    ungroup()


ggplot(binding_regions) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=4, max.overlaps = 100, nudge_y=-100, angle=90, force=5
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


binding_regions_label <- binding_regions %>% 
    filter(padj_peak<1e-4) %>% 
    group_by(start, end) %>% 
    top_n(3, estimate)

ggplot(binding_regions_label) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


plot_df <- binding_regions %>% 
    group_by(start, end) %>% 
    summarize(padj_peak=min(padj_peak))
p_mot <- ggplot(plot_df) +
    geom_text_repel(
        data = binding_regions_label,
        mapping = aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    ggnewscale::new_scale_color() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    xlim(c(start(region), end(region))) + no_legend()

p_cov <- CoveragePlot(
    atac_cov, 
    region = region, 
    group.by = 'lineage', 
    annotation = FALSE,
    peaks = FALSE
) & scale_fill_manual(values=lineage_colors) &
    theme_void() & no_legend() & theme(strip.text = element_blank())

p_an <- AnnotationPlot(
    atac_cov,
    region = region
) + theme_void() 

p_pe <- PeakPlot(
    atac_cov,
    region = region
) + theme_void() + no_legend()

p_cov / p_an / p_pe / p_mot + plot_layout(heights=c(5,1,1,3)) & theme(text = element_text(size=5))
ggsave('plots/paper/supp_grn/grn_reg_regions_neurod6_tracks.pdf', width=4.4, height=3.7, unit='cm')




#### DLX5 ####

to_tf <- 'DLX5'
annot <- gene_annot[gene_annot$gene_name==to_tf]
strand <- mode(as.character(strand(annot)))
# strand <- '+'
region <- GRanges(
    seqnames = as.character(seqnames(annot))[[1]], 
    ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
    strand = strand
)
# region <- resize(region, width=0, fix='start')
region <- Extend(region, upstream = 4000, downstream=10)


binding_regions <- modules %>% filter(target==to_tf) %>% pull(peaks) %>% str_split(';') %>% unlist() %>% StringToGRanges(sep=c('_', '_'))

binding_regions <- modules %>% filter(target==to_tf) %>% 
    separate_rows(peaks, sep=';') %>% 
    separate(peaks, sep='_', into=c('seqnames', 'start', 'end'), convert=T) %>% 
    filter(end>start(region), start<end(region)) %>% 
    group_by(start, end) %>% 
    mutate(padj_peak=min(padj)) %>% 
    mutate(padj_peak=pmax(padj, 1e-10)) %>% 
    mutate(pos_ratio=sum(sign(estimate)==1)/n()) %>% 
    ungroup()


ggplot(binding_regions) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=4, max.overlaps = 100, nudge_y=-100, angle=90, force=5
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


binding_regions_label <- binding_regions %>% 
    # filter(padj_peak<1e-4) %>% 
    group_by(start, end) %>% 
    top_n(3, -padj_peak)

ggplot(binding_regions_label) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


plot_df <- binding_regions %>% 
    group_by(start, end) %>% 
    summarize(padj_peak=min(padj_peak))
p_mot <- ggplot(plot_df) +
    geom_text_repel(
        data = binding_regions_label,
        mapping = aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    ggnewscale::new_scale_color() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    xlim(c(start(region), end(region))) + no_legend()

p_cov <- CoveragePlot(
    atac_cov, 
    region = region, 
    group.by = 'lineage', 
    annotation = FALSE,
    peaks = FALSE
) & scale_fill_manual(values=lineage_colors) &
    theme_void() & no_legend() & theme(strip.text = element_blank())

p_an <- AnnotationPlot(
    atac_cov,
    region = region
) + theme_void() 

p_pe <- PeakPlot(
    atac_cov,
    region = region
) + theme_void() + no_legend()

p_cov / p_an / p_pe / p_mot + plot_layout(heights=c(5,1,1,3)) & theme(text = element_text(size=5))
ggsave('plots/paper/supp_grn/grn_reg_regions_dlx6_tracks.pdf', width=4.4, height=3.7, unit='cm')




#### DLX5 ####

to_tf <- 'PAX6'
annot <- gene_annot[gene_annot$gene_name==to_tf]
strand <- mode(as.character(strand(annot)))
# strand <- '+'
region <- GRanges(
    seqnames = as.character(seqnames(annot))[[1]], 
    ranges = IRanges(start = min(start(annot)), end = max(end(annot))),
    strand = strand
)
# region <- resize(region, width=0, fix='start')
region <- Extend(region, upstream = 80000, downstream=0)

binding_regions <- modules %>% filter(target==to_tf) %>% pull(peaks) %>% str_split(';') %>% unlist() %>% StringToGRanges(sep=c('_', '_'))

binding_regions <- modules %>% filter(target==to_tf) %>% 
    separate_rows(peaks, sep=';') %>% 
    separate(peaks, sep='_', into=c('seqnames', 'start', 'end'), convert=T) %>% 
    filter(end>start(region), start<end(region)) %>% 
    group_by(start, end) %>% 
    mutate(padj_peak=min(padj)) %>% 
    mutate(padj_peak=pmax(padj, 1e-10)) %>% 
    mutate(pos_ratio=sum(sign(estimate)==1)/n()) %>% 
    ungroup()


ggplot(binding_regions) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=4, max.overlaps = 100, nudge_y=-100, angle=90, force=5
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


binding_regions_label <- binding_regions %>% 
    # filter(tf=='GLI3') %>%
    group_by(start, end) %>% 
    top_n(2, -padj_peak)

ggplot(binding_regions_label) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=10) +
    theme_void() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    ggnewscale::new_scale_color() +
    geom_text_repel(
        aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
        
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    xlim(c(start(region), end(region))) + no_legend()


plot_df <- binding_regions %>% 
    group_by(start, end) %>% 
    summarize(padj_peak=min(padj_peak))
p_mot <- ggplot(plot_df) +
    geom_text_repel(
        data = binding_regions_label,
        mapping = aes(x=(start+end)/2, y=0, color=factor(sign(estimate)), label=tf), 
        size=5/ggplot2::.pt, max.overlaps = 100, nudge_y=-100, angle=90, force=1,
    ) +
    scale_color_manual(values=c('#20877B', '#882160')) +
    ggnewscale::new_scale_color() +
    scale_color_gradientn(colors=greys(), limits=c(0, 10)) +
    geom_segment(aes(x=start, y=0, xend=end, yend=0, color=-log10(padj_peak)), size=4) +
    theme_void() +
    xlim(c(start(region), end(region))) + no_legend()

p_cov <- CoveragePlot(
    atac_cov, 
    region = region, 
    group.by = 'lineage', 
    annotation = FALSE,
    peaks = FALSE
) & scale_fill_manual(values=lineage_colors) &
    theme_void() & no_legend() & theme(strip.text = element_blank())

p_an <- AnnotationPlot(
    atac_cov,
    region = region
) + theme_void() 

p_pe <- PeakPlot(
    atac_cov,
    region = region
) + theme_void() + no_legend()

p_cov / p_an / p_pe / p_mot + plot_layout(heights=c(5,1,1,3)) & theme(text = element_text(size=5))
ggsave('plots/paper/supp_grn/grn_reg_regions_pax6_tracks.pdf', width=4.4, height=3.7, unit='cm')












