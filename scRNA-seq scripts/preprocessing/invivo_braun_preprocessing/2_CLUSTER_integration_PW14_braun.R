if(TRUE) {
  rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
  
  library(Seurat)
  dataloc <- "./"
  
  # Setup ----------------------------
  OutDir <- OutDirOrig <- dataloc
  setwd(OutDirOrig)

  invivo.obj <- readRDS("./braun.obj.processed.14.RDS")
  invivo.obj$sample <- "invivo"

  invitro.obj <- readRDS("./combined.obj.RDS")
  invitro.obj$sample <- "invitro"

  ls.Seurat <- c(invitro.obj, 
                 invivo.obj)
  
  features <- SelectIntegrationFeatures(object.list = ls.Seurat,
                                        nfeatures = 3000)
  
  RPS.genes <- !grepl("^RPL|^RS|^RPS|^MRP", features)
  MT.genes <- !grepl("^MT\\.|^MT-", features)
  MALAT1.genes <- !grepl("MALAT1", features)
  
  features <- features[RPS.genes & MT.genes & MALAT1.genes]
  
  anchors <- FindIntegrationAnchors(object.list = ls.Seurat,
                                    anchor.features = features)
  
  combined.obj <- IntegrateData(anchorset = anchors)
  
  DefaultAssay(combined.obj) <- "integrated"
  
  combined.obj <- CellCycleScoring(combined.obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  combined.obj <- ScaleData(combined.obj, vars.to.regress = c("S.Score", "G2M.Score"))
  combined.obj <- RunPCA(combined.obj)
  combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:50, seed.use = 1) 

  saveRDS(combined.obj, "integrated_braun_pw14.RDS")
} 
