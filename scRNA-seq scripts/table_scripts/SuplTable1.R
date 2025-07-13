rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)

# Setup ----------------------------
setwd("./")
DataDir <- "./scRNAseq_analysis/data/"

source("./helper_functions/colors_metadata.R")

combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))

Idents(combined.obj) <- combined.obj$annotation_fine_num
allmarkers <- FindAllMarkers(combined.obj,only.pos = TRUE)

saveRDS(allmarkers, paste0(DataDir,"allmarkers_RNA_snn_res.3.5.RDS"))

allmarkers_sub <- allmarkers[allmarkers$p_val_adj < 0.01 & 
                               allmarkers$avg_log2FC > 1,]

writexl::write_xlsx(allmarkers_sub, "tables/SuplTable1_DEG_scRNA_invitro.xlsx")
