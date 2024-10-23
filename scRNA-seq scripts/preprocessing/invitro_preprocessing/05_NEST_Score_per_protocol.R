rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(NESTscore)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
OutDir <- OutDirOrig <- "./scRNAseq_analysis/allinone/"
setwd(OutDirOrig)

DataDir <- "./scRNAseq_analysis/data/"
combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))

umap <- combined.obj@reductions$umap@cell.embeddings
df <- cbind(umap, combined.obj@meta.data)

pcadim <- 40
k0 <- 100
ratio_permissivness <- 0.75

Idents(combined.obj) <- combined.obj$protocol
ls.protocol <- SplitObject(combined.obj)

meta <- data.frame()

cell.line_cat_df <- matrix("cell line-driven", 
                           nrow = length(levels(combined.obj$cell.line)),
                           ncol = length(levels(combined.obj$protocol)))
rownames(cell.line_cat_df) <- levels(combined.obj$cell.line)
colnames(cell.line_cat_df) <- levels(combined.obj$protocol)

combined.obj@misc$protocol_driven <- rep(NA,4)
names(combined.obj@misc$protocol_driven) <- names(protocol_colors)[1:4]

thresh_list <- rep(NA,4)
names(thresh_list) <- levels(combined.obj$protocol)

for(obj in ls.protocol) {
  
  prot <- as.character(unique(obj$protocol))
  setwd(prot)
  
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:pcadim)
  
  NESTres <- NESTscore(obj, "cell.line", 
                       k_nn = k0, 
                       ndims = pcadim)
  obj <- NESTres$seuratobj
  thresh <- NEST_threshold_proposal(obj, "cell.line", 
                                    min_nsamples = 2)
  thresh_list[[prot]] <- thresh
  
  obj$protocol_driven <- c("cell line-driven","protocol-driven")[(obj@meta.data[,"NESTscore_cell.line"] > thresh) + 1]
  obj$protocol_driven <- factor(obj$protocol_driven, levels = c("cell line-driven","protocol-driven"))

  for(cl in unique(obj$cell.line)) { 
    ratios <- table(obj$protocol_driven[which(obj$cell.line == cl)])/
      sum(obj$cell.line == cl)
    if(ratios["protocol-driven"] >= ratio_permissivness) { 
      cell.line_cat_df[cl,prot] <- "protocol-driven"
    }
  }
  
  labels <- obj@meta.data %>% 
    group_by(annotation_coarse) %>% 
    summarise(proportion = round(sum(protocol_driven == "protocol-driven")/n(),2))

  meta <- rbind(meta, cbind("cells" = colnames(obj), 
                            "NESTscore_cell.line" = obj$NESTscore_cell.line, 
                            "protocol_driven" = obj$protocol_driven))
  
  setwd(OutDirOrig)
  
}

rownames(meta) <- meta$cells
meta$cells <- NULL
combined.obj <- AddMetaData(combined.obj, metadata = meta)
combined.obj$protocol_driven <- c("cell line-driven","protocol-driven")[as.numeric(combined.obj$protocol_driven)]

labels <- combined.obj@meta.data %>% 
  group_by(annotation_coarse) %>% 
  summarise(proportion = round(sum(protocol_driven == "protocol-driven")/n(),2))

write.csv(cell.line_cat_df, paste0(DataDir,"cell.line_overview_permissivness.csv"))

combined.obj@misc$protocol_driven <- thresh_list
saveRDS(combined.obj, paste0(DataDir,"combined.obj.RDS"))
