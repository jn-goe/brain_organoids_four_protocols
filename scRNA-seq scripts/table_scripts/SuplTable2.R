rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(dplyr)

# Setup ----------------------------
setwd("./")
DataDir <- "./scRNAseq_analysis/data/"

### in vitro ###
combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))

marker <- readRDS(paste0(DataDir,"allmarkers_RNA_snn_res.3.5.RDS"))
marker %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

conversion_table <- combined.obj@meta.data[, c("RNA_snn_res.3.5",
                                               "annotation_fine_num", 
                                               "annotation_fine", 
                                               "annotation_coarse")]

conversion_table <- conversion_table[!duplicated(conversion_table),]
rownames(conversion_table) <- conversion_table$annotation_fine_num
conversion_table$top10_marker <- ""

for(ct in conversion_table$annotation_fine_num) { 
  marker_ct <- top10$gene[top10$cluster == ct]
  conversion_table[ct,5] <- paste(marker_ct,collapse=", ")
}

rownames(conversion_table) <- NULL
rm(marker)

conversion_table_invitro <- conversion_table

### in vivo - Braun et al ###
combined.obj <- readRDS(paste0(DataDir,"integrated_braun_pw14.RDS"))
obj <- combined.obj[,which(combined.obj$sample == "invivo")]

conversion_table <- obj@meta.data[, c("Subregion",
                                      "CellClass",
                                      "annotation")]
conversion_table <- conversion_table[!duplicated(conversion_table),]
rownames(conversion_table) <- NULL
colnames(conversion_table)[1:2] <- paste0(colnames(conversion_table)[1:2], " (original)")
colnames(conversion_table)[3] <- paste0("combined ",colnames(conversion_table)[3], " used for visualizations")
conversion_table_braun <- conversion_table

### in vivo - Bhaduri et al ###
combined.obj <- readRDS(paste0(DataDir,"integrated_bhaduri_GW141718_sampled.RDS"))
obj <- combined.obj[,which(combined.obj$sample == "invivo")]

conversion_table <- obj@meta.data[, c("structure",
                                               "cell_type",
                                               "annotation")]
conversion_table <- conversion_table[!duplicated(conversion_table),]
rownames(conversion_table) <- NULL
colnames(conversion_table)[1:2] <- paste0(colnames(conversion_table)[1:2], " (original)")
colnames(conversion_table)[3] <- paste0("combined ",colnames(conversion_table)[3], " used for visualizations")
conversion_table_bhaduri <- conversion_table

conversion_table_list <- list("invitro" = conversion_table_invitro,
                              "braun_etal" = conversion_table_braun,
                              "bhaduri_etal" = conversion_table_bhaduri)

writexl::write_xlsx(conversion_table_list, 
                    "tables/SuplTable2_clustering_hierarchy.xlsx")
