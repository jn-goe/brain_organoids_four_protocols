rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)

# Setup ----------------------------
setwd("./")
DataDir <- "./scRNAseq_analysis/data/"

combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))

Idents(combined.obj) <- combined.obj$annotation_fine_num
allmarkers <- FindAllMarkers(combined.obj,only.pos = TRUE)

writexl::write_xlsx(allmarkers, "tables/SuplTable1_DEG_scRNA_invitro.xlsx")
saveRDS(allmarkers, paste0(DataDir,"allmarkers_RNA_snn_res.3.5.RDS"))
