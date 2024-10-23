rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(dplyr)
library(ggplot2)
library(Seurat)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./bulkRNAseq_analysis/data/"
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

### Sample Overview ###
percellineandprotocol <- metadata %>% group_by(protocol, cell.line, 
                                              repetition, depth, time_str, time_num) %>%  tally()
percellineandprotocol$n <- NULL
percellineandprotocol$replicate <- paste0("ORG",percellineandprotocol$repetition)
percellineandprotocol$repetition <- NULL

percellineandprotocol$replicate <- factor(percellineandprotocol$replicate, levels = c("ORG1","ORG2","ORG3"))

percellineandprotocol$cell.line <- factor(percellineandprotocol$cell.line, 
                                          levels = names(cell.line_colors))
percellineandprotocol$cell.line_num <- as.numeric(percellineandprotocol$cell.line) + 1/(length(unique(as.numeric(percellineandprotocol$replicate)))+1)*as.numeric(percellineandprotocol$replicate)

percellineandprotocol$time_str <- factor(percellineandprotocol$time_str, 
                                         levels = sort(unique(percellineandprotocol$time_str)))
percellineandprotocol$time <- as.numeric(percellineandprotocol$time_str) + 1/(length(unique(percellineandprotocol$protocol))+1)*as.numeric(percellineandprotocol$protocol)

p1 <- ggplot() + 
  geom_point(data = percellineandprotocol[percellineandprotocol$depth > 100000,],
             aes(y = time, 
                 x = cell.line_num, 
                 col = protocol)) + 
  scale_color_manual(values=protocol_colors) + xlab("") + ylab("") + 
  theme_minimal() + 
  scale_x_continuous(limits = range(percellineandprotocol$cell.line_num),
                     breaks = 1:length(levels(percellineandprotocol$cell.line)) + .5,
                     labels = levels(percellineandprotocol$cell.line)) +
  scale_y_continuous(limits = range(percellineandprotocol$time),
                     breaks = 1:length(levels(percellineandprotocol$time_str)) + .5,
                     labels = levels(percellineandprotocol$time_str)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + NoLegend()
p1

try(p1 + ggsave("figures/SFig4/SuplFig4A_sample_overview.pdf", width = 3, height = 3), silent = T)

