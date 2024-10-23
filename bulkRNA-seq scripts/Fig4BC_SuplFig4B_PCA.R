rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(dplyr)
library(ggplot2)
library(Seurat)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./bulkRNAseq_analysis/data/"
dir.create("figures/Fig4", showWarnings = FALSE)
dir.create("figures/SFig4", showWarnings = FALSE)

cell.line_colors <- cell.line_colors[-1]

#### DATA LOADING AND PREPROCESSING
gene_symbols = read.table(paste0(DataDir,"ii.GRCh38_202211_mapTable.tsv"), header=TRUE)
colnames(gene_symbols) = c('ensembl', 'symbol')

metadata <- read.table(paste0(DataDir,"bulkRNAseq_metadata.csv"), sep=',', header=TRUE)
expr <- readRDS(paste0(DataDir,"allsamples_counts_TPM.RDS"))
counts <- readRDS(paste0(DataDir,"allsamples_counts.RDS"))

# EXCLUDE SAMPLES WITH LESS THAN 100,000 READS
metadata$depth = apply(counts, 2, sum)
expr = expr[,metadata$depth > 100000]
counts = counts[,metadata$depth > 100000]
metadata = metadata[metadata$depth > 100000,]

metadata$cell.line <- factor(metadata$cell.line, levels = names(cell.line_colors))
metadata$protocol <- factor(metadata$protocol, levels = names(protocol_colors))

######## PCA ######## 
rsum = apply(expr,1,sum)
expr = expr[rsum >= 10,]
expr = log2(expr+1)

vars = apply(expr, 1, var)
ind = order(vars, decreasing=TRUE)[1:2000]

pc = prcomp(t(expr[ind,]), scale.=TRUE, center=TRUE)
pc$expl_var = (pc$sdev / sum(pc$sdev)) * 100

for(i in seq(1,8,2)){
  metadata$pc_a = pc$x[,i]
  metadata$pc_b = pc$x[,i+1]
  
  if(i == 1) {
    p_time <- (ggplot(data=metadata, aes(x=pc_a, y=pc_b, color=time_num)) +
                 geom_point(size=5) +
                 xlab(paste0('PC ', i, ' (Expl. Var. = ', round(pc$expl_var[i], 1), ' %)')) + 
                 ylab(paste0('PC ', i+1, ' (Expl. Var. = ', round(pc$expl_var[i+1], 1), ' %)')) +
                 theme_minimal() +
                 theme(text=element_text(size=25)))
    try(p_time + ggsave(paste0("figures/Fig4/Fig4B_PCA12_time.pdf"), width = 9, height = 8),silent = T)
  } else if(i == 3) {
    p_prot <- (ggplot(data=metadata, aes(x=pc_a, y=pc_b, color=protocol)) +
                 geom_point(size=5) +
                 xlab(paste0('PC ', i, ' (Expl. Var. = ', round(pc$expl_var[i], 1), ' %)')) + 
                 ylab(paste0('PC ', i+1, ' (Expl. Var. = ', round(pc$expl_var[i+1], 1), ' %)')) +
                 theme_minimal() +
                 theme(text=element_text(size=25)) +
                 scale_color_manual(values=protocol_colors))
    try(p_prot + ggsave(paste0("figures/Fig4/Fig4C_PCA34_protocol.pdf"), width = 9, height = 8),silent = T)
  } else if(i == 7) {
    p_cell <- (ggplot(data=metadata, aes(x=pc_a, y=pc_b, color=cell.line)) +
                 geom_point(size=5) +
                 xlab(paste0('PC ', i, ' (Expl. Var. = ', round(pc$expl_var[i], 1), ' %)')) + 
                 ylab(paste0('PC ', i+1, ' (Expl. Var. = ', round(pc$expl_var[i+1], 1), ' %)')) +
                 theme_minimal() +
                 theme(text=element_text(size=25)) +
                 scale_color_manual(values=cell.line_colors))
    try(p_cell + ggsave(paste0("figures/SFig4/SuplFig4B_PCA78_cell.line.pdf"), width = 9, height = 8),silent = T)
  }
}
