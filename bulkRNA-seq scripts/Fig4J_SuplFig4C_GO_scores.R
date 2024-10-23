rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(dplyr)
library(ggplot2)
library(Seurat)

source("./helper_functions/colors_metadata.R")
source("./helper_functions/bulk_functions.R")

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

rownames(expr) <- rownames(counts) <-
  make.unique(map_symbols(rownames(expr), gene_symbols))

bulk_obj <- CreateSeuratObject(counts, meta.data = metadata, min.cells = 10)
bulk_obj <- NormalizeData(bulk_obj) # log(expr/100+1)
bulk_obj <- FindVariableFeatures(bulk_obj, nfeatures = 2000)
bulk_obj <- ScaleData(bulk_obj)
bulk_obj <- RunPCA(bulk_obj)
bulk_obj <- RunUMAP(bulk_obj, dims = 1:50, n.components = 2L)

saveRDS(bulk_obj, paste0(DataDir,"bulk_obj.RDS"))

ensembl <- biomaRt::useEnsembl("ensembl",
                                dataset = "hsapiens_gene_ensembl",
                                mirror = NULL)

GOlist <- list(go1="GO:0006096", # Glycolysis
               go2="GO:0034976" # ER-stress
               )

for(GO in GOlist) { 
  obj <- bulk_obj
  
  suppressMessages(suppressWarnings(genes <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                                                   GO, c("SYMBOL"),"GOALL")$SYMBOL %>% unique()))
  
  GO_name <- Term(GOTERM)[which(names(Term(GOTERM)) == GO)]
  
  genes <- intersect(rownames(obj), genes)
  obj <- Seurat::AddModuleScore(obj, features = list(genes), name = "go_score")
  
  range_goterm = range(obj$go_score1)
  
  ##### dorsal #####
  
  mat = matrix(NA,
               nrow = length(unique(obj$time_str)),
               ncol = length(levels(obj$cell.line)))
  colnames(mat) = levels(obj$cell.line)
  rownames(mat) = unique(obj$time_str)
  
  for(t in unique(obj$time_str)) {
    
    if(t == "day000") {
      df = obj@meta.data[obj$time_str == t,] %>%
        group_by(cell.line) %>%
        summarise(across(go_score1, mean, na.rm = TRUE))
      
      mat[as.character(t),as.character(df$cell.line)] = df$go_score1
    } else {
      df = obj@meta.data[obj$time_str == t & obj$protocol == "dorsal",] %>%
        group_by(cell.line) %>%
        summarise(across(go_score1, mean, na.rm = TRUE))
      
      mat[as.character(t),as.character(df$cell.line)] = df$go_score1
    }
  }
  
  mat = mat[c("day000","day013-016","day025","day040","day080","day120"),]
  rownames(mat)[grepl("day013",rownames(mat))] = "day013/16"
  
  pdf(paste0("figures/Fig4/Fig4J_",sub(":","",GO),"_GOscore_overtime_cellline_dorsal.pdf"),
      width = 3.6, height = 3)
  pheatmap::pheatmap(t(mat),
                     cluster_rows = F,
                     cluster_cols = F,
                     main =paste0("dorsal protocol"),
                     color = colorRampPalette(c("darkblue","grey70","darkred"))(20), na_col = "transparent",
                     breaks = round(seq(range_goterm[1], range_goterm[2], length.out = 20),3))
  dev.off()  
  
  ##### ventral #####
  
  mat = matrix(NA,
               nrow = length(unique(obj$time_str)),
               ncol = length(levels(obj$cell.line)))
  colnames(mat) = levels(obj$cell.line)
  rownames(mat) = unique(obj$time_str)
  
  for(t in unique(obj$time_str)) {
    
    if(t == "day000") {
      df = obj@meta.data[obj$time_str == t,] %>%
        group_by(cell.line) %>%
        summarise(across(go_score1, mean, na.rm = TRUE))
      
      mat[as.character(t),as.character(df$cell.line)] = df$go_score1
    } else {
      df = obj@meta.data[obj$time_str == t & obj$protocol == "ventral",] %>%
        group_by(cell.line) %>%
        summarise(across(go_score1, mean, na.rm = TRUE))
      
      mat[as.character(t),as.character(df$cell.line)] = df$go_score1
    }
  }
  
  mat = mat[c("day000","day013-016","day025","day040","day080","day120"),]
  rownames(mat)[grepl("day013",rownames(mat))] = "day013/16"
  
  pdf(paste0("figures/SFig4/SuplFig4C_",sub(":","",GO),"_GOscore_overtime_cellline_ventral.pdf"),
      width = 3.6, height = 3)
  pheatmap::pheatmap(t(mat),
                     cluster_rows = F,
                     cluster_cols = F,
                     main =paste0("ventral protocol"),
                     color = colorRampPalette(c("darkblue","grey70","darkred"))(20), na_col = "transparent",
                     breaks = round(seq(range_goterm[1], range_goterm[2], length.out = 20),3))
  dev.off()  
  ##### midbrain #####
  
  mat = matrix(NA,
               nrow = length(unique(obj$time_str)),
               ncol = length(levels(obj$cell.line)))
  colnames(mat) = levels(obj$cell.line)
  rownames(mat) = unique(obj$time_str)
  
  for(t in unique(obj$time_str)) {
    
    if(t == "day000") {
      df = obj@meta.data[obj$time_str == t,] %>%
        group_by(cell.line) %>%
        summarise(across(go_score1, mean, na.rm = TRUE))
      
      mat[as.character(t),as.character(df$cell.line)] = df$go_score1
    } else {
      df = obj@meta.data[obj$time_str == t & obj$protocol == "midbrain",] %>%
        group_by(cell.line) %>%
        summarise(across(go_score1, mean, na.rm = TRUE))
      
      mat[as.character(t),as.character(df$cell.line)] = df$go_score1
    }
  }
  
  mat = mat[c("day000","day013-016","day025","day040","day080","day120"),]
  rownames(mat)[grepl("day013",rownames(mat))] = "day013/16"
  
  pdf(paste0("figures/SFig4/SuplFig4C_",sub(":","",GO),"_GOscore_overtime_cellline_midbrain.pdf"),
      width = 3.6, height = 3)
  pheatmap::pheatmap(t(mat),
                     cluster_rows = F,
                     cluster_cols = F,
                     main =paste0("midbrain protocol"),
                     color = colorRampPalette(c("darkblue","grey70","darkred"))(20), na_col = "transparent",
                     breaks = round(seq(range_goterm[1], range_goterm[2], length.out = 20),3))
  dev.off()
  
  
  ##### striatum #####
  
  mat = matrix(NA,
               nrow = length(unique(obj$time_str)),
               ncol = length(levels(obj$cell.line)))
  colnames(mat) = levels(obj$cell.line)
  rownames(mat) = unique(obj$time_str)
  
  for(t in unique(obj$time_str)) {
    
    if(t == "day000") {
      df = obj@meta.data[obj$time_str == t,] %>%
        group_by(cell.line) %>%
        summarise(across(go_score1, mean, na.rm = TRUE))
      
      mat[as.character(t),as.character(df$cell.line)] = df$go_score1
    } else {
      df = obj@meta.data[obj$time_str == t & obj$protocol == "striatum",] %>%
        group_by(cell.line) %>%
        summarise(across(go_score1, mean, na.rm = TRUE))
      
      mat[as.character(t),as.character(df$cell.line)] = df$go_score1
    }
  }
  
  mat = mat[c("day000","day013-016","day025","day040","day080","day120"),]
  rownames(mat)[grepl("day013",rownames(mat))] = "day013/16"
  
  pdf(paste0("figures/SFig4/SuplFig4C_",sub(":","",GO),"_GOscore_overtime_cellline_striatum.pdf"),
      width = 3.6, height = 3)
  pheatmap::pheatmap(t(mat),
                     cluster_rows = F,
                     cluster_cols = F,
                     main = paste0("striatum protocol"),
                     color = colorRampPalette(c("darkblue","grey70","darkred"))(20), na_col = "transparent",
                     breaks = round(seq(range_goterm[1], range_goterm[2], length.out = 20),3))
  dev.off()
  
}
