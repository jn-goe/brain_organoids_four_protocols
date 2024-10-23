# Setup ----------------------------

setwd('./bulkRNAseq_analysis')

DataDirsc <- "./scRNAseq_analysis/data/"
DataDirbulk <- "./bulkRNAseq_analysis/data/"

metadata <- read.csv(paste0(DataDirbulk,"bulkRNAseq_metadata.csv"))

cell.line_cat_df <- read.csv(paste0(DataDirsc,"cell.line_overview_permissivness.csv"))
rownames(cell.line_cat_df) <- cell.line_cat_df$X
cell.line_cat_df$X <- NULL

cell.line_driven <- paste0(rownames(cell.line_cat_df)[which(cell.line_cat_df == "cell line-driven", arr.ind = T)[,1]],"-",
                           colnames(cell.line_cat_df)[which(cell.line_cat_df == "cell line-driven", arr.ind = T)[,2]])

metadata$growth <- "protocol-driven"
metadata$growth[which(paste0(metadata$cell.line,"-",metadata$protocol) %in% cell.line_driven)] <- "cell line-driven"
  
write.csv(metadata, paste0(DataDirbulk,"bulkRNAseq_metadata.csv"))
