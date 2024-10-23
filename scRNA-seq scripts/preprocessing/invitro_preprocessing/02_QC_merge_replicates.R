rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(dplyr)
library(reticulate)
scrublet <- import("scrublet")

source("./helper_functions/functions_QC.R")

# Setup ----------------------------
OutDirOrig <- "./scRNAseq_analysis/separate_organoids/"
setwd(OutDirOrig)

folderlist <- list.files(paste0(OutDirOrig,"cellranger"))
folderlist <- folderlist[grepl(".Rds", folderlist)]
org <- list(org1=folderlist[grepl("1",substr(folderlist,3,3))],
            org2=folderlist[grepl("2",substr(folderlist,3,3))])

countorg2 <- 1

colnames <- c("organoid", "ncells_beforeQC",	"ncells_afterQC", "median.UMI.counts_afterQC",	"median.gene.counts_afterQC")
QC_table <- matrix(nrow = length(org$org1) + length(org$org1), ncol = length(colnames))
colnames(QC_table) <- colnames
iteration <- 1

mito.range_for_all = c(0.001, 0.15) 
ribo.range_for_all = c(0.001, 0.5)
nfeature.range_for_protocol = c(500, 5000)
ncount.range_for_protocol = c(500, 10000)

max_ratio_gene = 0.1

for(countorg1 in 1:length(org$org1)) {
  
  protocol <- substr(org$org1[[countorg1]],2,2)

  if(substr(org$org1[[countorg1]],1,2) == substr(org$org2[[countorg2]],1,2)) {
    
    obj.org1 <- readRDS(paste0(OutDirOrig,"cellranger/",org$org1[[countorg1]]))
    obj.org1 <- RenameCells(obj.org1, add.cell.id = unique(obj.org1$orig.ident))
    
    cellnumber_beforeQC1 <- dim(obj.org1)[2]
    
    obj.org2 <- readRDS(paste0(OutDirOrig,"cellranger/",org$org2[[countorg2]]))
    obj.org2 <- RenameCells(obj.org2, add.cell.id = unique(obj.org2$orig.ident))
    
    cellnumber_beforeQC2 <- dim(obj.org2)[2]
    
    setwd(paste0("organoids_",substr(org$org1[[countorg1]],1,2)))
  
    #-------------------------------------#
    
    obj.org1$ratio_max_gene <- as.matrix(apply(obj.org1@assays$RNA@layers$counts, 2, fun_max_ratio))
    obj.org1$ratio_max_gene_over_thresh <- obj.org1$ratio_max_gene > max_ratio_gene
    Idents(obj.org1) <- obj.org1$ratio_max_gene_over_thresh
    obj.org1 <- subset(obj.org1, ident = F)
    
    new_counts <- obj.org1@assays$RNA@layers$counts[which(rownames(obj.org1)!="MALAT1"),]
    rownames(new_counts) <- rownames(obj.org1)[which(rownames(obj.org1)!="MALAT1")]
    colnames(new_counts) <- colnames(obj.org1)
    
    obj.org1 <- CreateSeuratObject(new_counts, meta.data = obj.org1@meta.data)
    obj.org1 <- NormalizeData(obj.org1)
    
    obj.org1[["percent.mito"]] <- PercentageFeatureSet(obj.org1, pattern = "^MT\\.|^MT-")
    obj.org1[["percent.mito"]] <- obj.org1[["percent.mito"]]/100
    
    obj.org1[["percent.ribo"]] <- PercentageFeatureSet(obj.org1, pattern = "^RPL|^RS|^RPS|^MRP")
    obj.org1[["percent.ribo"]] <- obj.org1[["percent.ribo"]]/100
    
    QC_res <- plot_QC(obj.org1,
                      nfeature.range = nfeature.range_for_protocol,
                      ncount.range = ncount.range_for_protocol,
                      mito.range = mito.range_for_all,
                      ribo.range = ribo.range_for_all)

    obj.org1 <- QC_res[[2]]
    Idents(obj.org1) <- obj.org1$passQC
    
    scrub = scrublet$Scrublet(t(as.matrix(obj.org1@assays$RNA@layers$counts)))
    scrub_doublets = scrub$scrub_doublets()
    predicted_doublets = scrub_doublets[[2]] 
    names(predicted_doublets) <- colnames(obj.org1)
    
    cells.no.doublets <- names(predicted_doublets[!predicted_doublets])
    obj.org1 <- subset(obj.org1, cells = cells.no.doublets)
    
    cells_pass_QC <- automatic_QC(obj.org1, 
                                  nFeature_RNA.range = nfeature.range_for_protocol,
                                  nCount_RNA.range = ncount.range_for_protocol,
                                  mito.range = mito.range_for_all,
                                  ribo.range = ribo.range_for_all)
    
    obj.org1 <- subset(obj.org1, cells = cells_pass_QC)
    
    QC_table[iteration,] <- c(unique(obj.org1$orig.ident),
                              cellnumber_beforeQC1,
                              sum(obj.org1$replicate == "ORG1"),
                              round(median(obj.org1[,which(obj.org1$replicate == "ORG1")]$nCount_RNA),1),
                              round(median(obj.org1[,which(obj.org1$replicate == "ORG1")]$nFeature_RNA),1))
    
    #-------------------------------------#
    
    obj.org2$ratio_max_gene <- as.matrix(apply(obj.org2@assays$RNA@layers$counts, 2, fun_max_ratio))
    obj.org2$ratio_max_gene_over_thresh <- obj.org2$ratio_max_gene > max_ratio_gene
    Idents(obj.org2) <- obj.org2$ratio_max_gene_over_thresh
    obj.org2 <- subset(obj.org2, ident = F)
    
    new_counts <- obj.org2@assays$RNA@layers$counts[which(rownames(obj.org2)!="MALAT1"),]
    rownames(new_counts) <- rownames(obj.org2)[which(rownames(obj.org2)!="MALAT1")]
    colnames(new_counts) <- colnames(obj.org2)
    
    obj.org2 <- CreateSeuratObject(new_counts, meta.data = obj.org2@meta.data)
    obj.org2 <- NormalizeData(obj.org2)
    
    obj.org2[["percent.mito"]] <- PercentageFeatureSet(obj.org2, pattern = "^MT\\.|^MT-")
    obj.org2[["percent.mito"]] <- obj.org2[["percent.mito"]]/100
    
    obj.org2[["percent.ribo"]] <- PercentageFeatureSet(obj.org2, pattern = "^RPL|^RS|^RPS|^MRP")
    obj.org2[["percent.ribo"]] <- obj.org2[["percent.ribo"]]/100
    
    QC_res <- plot_QC(obj.org2,
                      nfeature.range = nfeature.range_for_protocol,
                      ncount.range = ncount.range_for_protocol,
                      mito.range = mito.range_for_all,
                      ribo.range = ribo.range_for_all)

    obj.org2 <- QC_res[[2]]
    Idents(obj.org2) <- obj.org2$passQC
    
    scrub = scrublet$Scrublet(t(as.matrix(obj.org2@assays$RNA@layers$counts)))
    scrub_doublets = scrub$scrub_doublets()
    predicted_doublets = scrub_doublets[[2]]
    names(predicted_doublets) <- colnames(obj.org2)
    
    cells.no.doublets <- names(predicted_doublets[!predicted_doublets])
    obj.org2 <- subset(obj.org2, cells = cells.no.doublets)
    
    cells_pass_QC <- automatic_QC(obj.org2, 
                                  nFeature_RNA.range = nfeature.range_for_protocol,
                                  nCount_RNA.range = ncount.range_for_protocol,
                                  mito.range = mito.range_for_all,
                                  ribo.range = ribo.range_for_all)
    
    obj.org2 <- subset(obj.org2, cells = cells_pass_QC)
    
    iteration <- iteration + 1
    
    QC_table[iteration,] <- c(unique(obj.org2$orig.ident),
                              cellnumber_beforeQC2,
                              sum(obj.org2$replicate == "ORG2"),
                              round(median(obj.org2[,which(obj.org2$replicate == "ORG2")]$nCount_RNA),1),
                              round(median(obj.org2[,which(obj.org2$replicate == "ORG2")]$nFeature_RNA),1))
    
    #-------------------------------------#
    # merge without batch correction
    combined.obj <- merge(x = obj.org1, y = obj.org2)
    combined.obj <- JoinLayers(combined.obj)
    
    iteration <- iteration + 1
    
    print(head(colnames(combined.obj)))
    
    saveRDS(combined.obj, "combined.obj.RDS")
    
    setwd(OutDirOrig) 
    
    countorg2 <- countorg2 + 1
  } else {
    
    #-------------------------------------#
    
    obj.org1 <- readRDS(paste0(OutDirOrig,"cellranger/",org$org1[[countorg1]]))
    obj.org1 <- RenameCells(obj.org1, add.cell.id = unique(obj.org1$orig.ident))
    
    cellnumber_beforeQC1 <- dim(obj.org1)[2]
    
    setwd(paste0("organoids_",substr(org$org1[[countorg1]],1,2)))
    
    obj.org1$ratio_max_gene <- as.matrix(apply(obj.org1@assays$RNA@layers$counts, 2, fun_max_ratio))
    obj.org1$ratio_max_gene_over_thresh <- obj.org1$ratio_max_gene > max_ratio_gene
    Idents(obj.org1) <- obj.org1$ratio_max_gene_over_thresh
    obj.org1 <- subset(obj.org1, ident = F)
    
    new_counts <- obj.org1@assays$RNA@layers$counts[which(rownames(obj.org1)!="MALAT1"),]
    rownames(new_counts) <- rownames(obj.org1)[which(rownames(obj.org1)!="MALAT1")]
    colnames(new_counts) <- colnames(obj.org1)
    
    obj.org1 <- CreateSeuratObject(new_counts, meta.data = obj.org1@meta.data)
    obj.org1 <- NormalizeData(obj.org1)
    
    obj.org1[["percent.mito"]] <- PercentageFeatureSet(obj.org1, pattern = "^MT\\.|^MT-")
    obj.org1[["percent.mito"]] <- obj.org1[["percent.mito"]]/100
    
    obj.org1[["percent.ribo"]] <- PercentageFeatureSet(obj.org1, pattern = "^RPL|^RS|^RPS|^MRP")
    obj.org1[["percent.ribo"]] <- obj.org1[["percent.ribo"]]/100
    
    QC_res <- plot_QC(obj.org1,
                      nfeature.range = nfeature.range_for_protocol,
                      ncount.range = ncount.range_for_protocol,
                      mito.range = mito.range_for_all,
                      ribo.range = ribo.range_for_all)

    obj.org1 <- QC_res[[2]]
    Idents(obj.org1) <- obj.org1$passQC
    
    scrub = scrublet$Scrublet(t(as.matrix(obj.org1@assays$RNA@layers$counts)))
    scrub_doublets = scrub$scrub_doublets()
    predicted_doublets = scrub_doublets[[2]] 
    names(predicted_doublets) <- colnames(obj.org1)
    
    cells.no.doublets <- names(predicted_doublets[!predicted_doublets])
    obj.org1 <- subset(obj.org1, cells = cells.no.doublets)
    
    cells_pass_QC <- automatic_QC(obj.org1, 
                                  nFeature_RNA.range = nfeature.range_for_protocol,
                                  nCount_RNA.range = ncount.range_for_protocol,
                                  mito.range = mito.range_for_all,
                                  ribo.range = ribo.range_for_all)
    
    obj.org1 <- subset(obj.org1, cells = cells_pass_QC)
    
    iteration <- iteration + 1
    
    QC_table[iteration,] <- c(unique(obj.org1$orig.ident),
                              cellnumber_beforeQC1,
                              sum(obj.org1$replicate == "ORG1"),
                              round(median(obj.org1[,which(obj.org1$replicate == "ORG1")]$nCount_RNA),1),
                              round(median(obj.org1[,which(obj.org1$replicate == "ORG1")]$nFeature_RNA),1)
    )
    
    print(head(colnames(obj.org1)))
    
    saveRDS(obj.org1, "combined.obj.RDS")
    setwd(OutDirOrig)
  }
}
