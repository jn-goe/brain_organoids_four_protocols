rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)

# Setup ----------------------------
OutDir <- OutDirOrig <- "./scRNAseq_analysis/bhaduri/data/"
setwd(OutDir)

folders <- list.files(OutDir)
folders <- folders[!grepl(pattern = ".Rds", x = folders)]

for(path in folders) {
  mat <- Read10X(paste0(path,"/",path,"/"))
  HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(mat), 
                                               unmapped.as.na = FALSE, 
                                               map = NULL, 
                                               species = "human")
  rownames(mat) <- make.unique(HGNC.updated$Suggested.Symbol)
  obj <- CreateSeuratObject(mat)
  
  obj$orig.ident <- path
  obj$time <- unlist(lapply(path, function(x) strsplit(x, "_")[[1]][1]))

  if(unlist(lapply(path, function(x) strsplit(x, "_")[[1]][2])) == "2") {
    obj$replicate <- 2
    obj$area <- unlist(lapply(path, function(x) strsplit(x, "_")[[1]][3]))
  } else {
    obj$replicate <- 1
    obj$area <- unlist(lapply(path, function(x) strsplit(x, "_")[[1]][2]))
  }
  
  saveRDS(obj, paste0(path,".Rds"))
}
