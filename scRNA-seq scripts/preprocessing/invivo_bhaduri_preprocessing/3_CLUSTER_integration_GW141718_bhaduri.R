if(TRUE) {
  rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
  
  library(Seurat)
  dataloc <- "./"
  
  # Setup ----------------------------
  OutDir <- OutDirOrig <- dataloc
  setwd(OutDirOrig)
  
  metadata.bhaduri <- readxl::read_xlsx("41586_2021_3910_MOESM3_ESM.xlsx", sheet = 1)
  metadata.bhaduri$cell.name <- sub("gw", "GW", metadata.bhaduri$cell.name)
  metadata.bhaduri$cell.name <- paste0(metadata.bhaduri$cell.name,"-1")
  metadata.bhaduri <- as.data.frame(metadata.bhaduri)
  rownames(metadata.bhaduri) <- metadata.bhaduri$cell.name
  
  invivo.ls <- readRDS("./bhaduri_GW141718_sampled_object_list.RDS")
  
  invivo.obj <- merge(invivo.ls[[1]],invivo.ls[2:length(invivo.ls)])
  invivo.obj <- JoinLayers(invivo.obj)
  invivo.obj$sample <- "invivo"
  invivo.obj <- AddMetaData(invivo.obj, metadata = metadata.bhaduri[colnames(invivo.obj),])
   
  invivo.ls <- SplitObject(invivo.obj, split.by = "individual")
  
  invitro.obj <- readRDS("./combined.obj.RDS")
  invitro.obj$sample <- "invitro"
   
  ls.Seurat <- c(invitro.obj,
                 invivo.ls)

  features <- SelectIntegrationFeatures(object.list = ls.Seurat,
                                        nfeatures = 3000)
  
  RPS.genes <- !grepl("^RPL|^RS|^RPS|^MRP", features)
  MT.genes <- !grepl("^MT\\.|^MT-", features)
  MALAT1.genes <- !grepl("MALAT1", features)

  features <- features[RPS.genes & MT.genes & MALAT1.genes]
  
  anchors <- FindIntegrationAnchors(object.list = ls.Seurat,
                                    anchor.features = features, 
                                    reference = 1)

  combined.obj <- IntegrateData(anchorset = anchors)
  
  DefaultAssay(combined.obj) <- "integrated"
  
  combined.obj <- CellCycleScoring(combined.obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  combined.obj <- ScaleData(combined.obj, vars.to.regress = c("S.Score", "G2M.Score"))
  combined.obj <- RunPCA(combined.obj)
  combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:50) 
  
  saveRDS(combined.obj, "integrated_bhaduri_GW141718_sampled.RDS")
} 
