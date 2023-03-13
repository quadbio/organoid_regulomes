##########################################################
############# COLLECTION OF MARKER GENES #################
##########################################################

mesenchyme_markers <- c('DCN', 'COL5A1', 'LUM')
colloid_markers <- c('TTR')
retina_markers <- c('VSX2', 'SIX6')
offtarget_markers <- c(mesenchyme_markers, colloid_markers, retina_markers)

main_markers <- c('PAX6', 'GLI3', 'SOX2', 'GSX2', 'NKX2-1', 'DLX2', 'GAD1', 'SLC32A1',
                  'BCL11B', 'TBR1', 'NEUROD6', 'MKI67', 'FOXG1', 'TTR', 'LHX9', 'TUBB3',
                  'DCN', 'EMX1', 'RSPO3', 'MEIS2', 'HOXB2', 'OTX2', 'ZIC1', 'WLS', 'SIX3')

early_markers <- c('POU5F1', 'LIN28A', 'KRT8')
ventral_markers <- c('GSX2', 'DLX5', 'GAD1')
midhind_markers <- c('LHX1', 'LHX5', 'LHX9', 'RSPO3')
cortical_markers <- c('GLI3', 'EOMES', 'NEUROD6')
lineage_markers <- c(early_markers, ventral_markers, midhind_markers, cortical_markers)

all_markers <- union(union(main_markers, lineage_markers), offtarget_markers)

neuron_markers <- unique(c('MAP2', 'ENO2', 'FOXG1', 'NEUROD6', 'NEUROD2', 'SLC17A7',
    'DLX5', 'DLX2', 'DLX1', 'GAD1', 'GAD2', 'SLC32A1', 'ISL1', 'EBF1', 'NHLH2', 'SLC17A6',
    'LHX9', 'TFAP2B', 'GBX2', 'SHOX2', 'PCP4', 'RSPO3', 'RELN', 'GATA3', 'OTX2',
    'SOX14', 'GLI3', 'EOMES', 'PAX3', 'PAX8', 'LHX9', 'PCP4','OTX2', 'GATA3',
    'DLX1', 'SOX6', 'VIP', 'ISL1', 'ONECUT3'))

patt_markers <- unique(c('SHH', 'RSPO2', 'WNT1', 'WNT3A', 'WNT8B', 'WNT9A', 'CNPY1',
    'FGF17','FGF8', 'TRH', 'DLK1', 'LMX1A', 'NTN1', 'PHOX2', 'MSX3', 'SOST', 'TAL2',
    'FGF15', 'BMP7', 'BMP4', 'WNT9A', 'WNT1', 'FGF15', 'NRG3', 'NRG1',
    'IRX3', 'IRX5', 'SIX3', 'VSX1', 'ZIC1', 'DMRT3', 'FEZF1', 'MSX1',
    'NKX2-1', 'DCX', 'BARHL2', 'TAL2'))


organiser_markers <- c(
    # Cortical Hem
    'SHH', 'RSPO2', 'WNT1', 'WNT3B', 'WNT8B', 'WNT9B',
    # Isthmus
    'CNPY1', 'FGF17', 'FGF8', 'TRH',
    # Other
    'DLK1', 'LMX1A', 'NTN1', 'PHOX2B', 'MSX2', 'SOST',
    'TAL2', 'WNT8B', 'FGF15', 'FGF8', 'TRH'
)


# Cell cycle markers from Zhisong
cc_genes <- readLines(con='~/resources/gene_sets/G2M_genes.txt')

# Cell cycle markers from Tirosh et al, 2015
cc_genes_tirosh <- readLines(con='~/resources/gene_sets/regev_lab_cell_cycle_genes.txt')

s_genes <- cc_genes[1:43]
g2m_genes <- cc_genes[44:97]

cc_genes_all <- union(cc_genes, cc_genes_tirosh)

#### RETINA ####
cone_markers <- c('ARR3','GUCA1A','PDE6H','GUCA1C')
rod_markers <- c('NRL','GNGT1','GNAT1','PDE6G','SAG','RHO','PDE6B','NR2E3','ROM1')
prpc_markers <- c('PDC','AIPL1','MIR7-3HG','EYS','CRX','PRDM1','NEUROD1','OTX2','RAX2')
bc_markers <- c('ISL1','VSX2','GRM6','GRIK1','VSX1','IRX6')
gc_markers <- c('CRABP1', 'GFAP', 'GLUL')
mg_markers <- c('RLBP1','SLC1A3','NEFL','NEFM')

retina_markers <- c(cone_markers, rod_markers, prpc_markers, bc_markers, gc_markers, mg_markers)
