if(TRUE) { 
  rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
  
  library(Seurat)
  library(ggplot2)
  library(NESTscore)
  
  combined.obj <- readRDS("integrated_braun_pw14.RDS")
  
  pcadim <- 50
  k0 <- 100
  
  ############## DORSAL ################

  cells_dorsal <- c(colnames(combined.obj)[which(combined.obj$protocol == "dorsal")],
                    colnames(combined.obj)[which(combined.obj$sample == "invivo")])
  
  dorsal <- subset(combined.obj, cells = cells_dorsal)
  dorsal <- RunPCA(dorsal)
  
  NESTres <- NESTscore(dorsal, "sample", 
                       k_nn = k0, 
                       ndims = pcadim)
  dorsal <- NESTres$seuratobj
  thresh <- NEST_threshold_proposal(dorsal, "sample", min_nsamples = 1)
  
  dorsal$invivo_mixed <- c("not well-mixed","well-mixed")[(dorsal@meta.data[,"NESTscore_sample"] > thresh) + 1]
  dorsal$invivo_mixed <- factor(dorsal$invivo_mixed, levels = c("not well-mixed","well-mixed"))
  
  meta_dorsal <- dorsal@meta.data[,c("NESTscore_sample","invivo_mixed")]
  rm(dorsal)
  
  ############## VENTRAL ################

  cells_ventral <- c(colnames(combined.obj)[which(combined.obj$protocol == "ventral")],
                    colnames(combined.obj)[which(combined.obj$sample == "invivo")])
  
  ventral <- subset(combined.obj, cells = cells_ventral)
  ventral <- RunPCA(ventral)
  
  NESTres <- NESTscore(ventral, "sample", 
                       k_nn = k0, 
                       ndims = pcadim)
  ventral <- NESTres$seuratobj
  thresh <- NEST_threshold_proposal(ventral, "sample", min_nsamples = 1)
  
  ventral$invivo_mixed <- c("not well-mixed","well-mixed")[(ventral@meta.data[,"NESTscore_sample"] > thresh) + 1]
  ventral$invivo_mixed <- factor(ventral$invivo_mixed, levels = c("not well-mixed","well-mixed"))
  
  meta_ventral <- ventral@meta.data[,c("NESTscore_sample","invivo_mixed")]
  rm(ventral)
  
  ############## MIDBRAIN ################

  cells_midbrain <- c(colnames(combined.obj)[which(combined.obj$protocol == "midbrain")],
                    colnames(combined.obj)[which(combined.obj$sample == "invivo")])
  
  midbrain <- subset(combined.obj, cells = cells_midbrain)
  midbrain <- RunPCA(midbrain)
  
  NESTres <- NESTscore(midbrain, "sample", 
                       k_nn = k0, 
                       ndims = pcadim)
  midbrain <- NESTres$seuratobj
  thresh <- NEST_threshold_proposal(midbrain, "sample", min_nsamples = 1)
  
  midbrain$invivo_mixed <- c("not well-mixed","well-mixed")[(midbrain@meta.data[,"NESTscore_sample"] > thresh) + 1]
  midbrain$invivo_mixed <- factor(midbrain$invivo_mixed, levels = c("not well-mixed","well-mixed"))
  
  meta_midbrain <- midbrain@meta.data[,c("NESTscore_sample","invivo_mixed")]
  rm(midbrain)
  
  ############## STRIATUM ################

  cells_striatum <- c(colnames(combined.obj)[which(combined.obj$protocol == "striatum")],
                    colnames(combined.obj)[which(combined.obj$sample == "invivo")])
  
  striatum <- subset(combined.obj, cells = cells_striatum)
  striatum <- RunPCA(striatum)
  
  NESTres <- NESTscore(striatum, "sample", 
                       k_nn = k0, 
                       ndims = pcadim)
  striatum <- NESTres$seuratobj
  thresh <- NEST_threshold_proposal(striatum, "sample", min_nsamples = 1)
  
  striatum$invivo_mixed <- c("not well-mixed","well-mixed")[(striatum@meta.data[,"NESTscore_sample"] > thresh) + 1]
  striatum$invivo_mixed <- factor(striatum$invivo_mixed, levels = c("not well-mixed","well-mixed"))
  
  meta_striatum <- striatum@meta.data[,c("NESTscore_sample","invivo_mixed")]
  rm(striatum)
  
  ############## ALL TOGETHER ################

  NESTres <- NESTscore(combined.obj, "sample", 
                       k_nn = k0, 
                       ndims = pcadim)
  combined.obj <- NESTres$seuratobj
  thresh <- NEST_threshold_proposal(combined.obj, "sample", min_nsamples = 1)
  
  combined.obj$invivo_mixed <- c("not well-mixed","well-mixed")[(combined.obj@meta.data[,"NESTscore_sample"] > thresh) + 1]
  combined.obj$invivo_mixed <- factor(combined.obj$invivo_mixed, levels = c("not well-mixed","well-mixed"))

  # store all results in object
  combined.obj$NESTscore_sample_dorsal <- 
    combined.obj$NESTscore_sample_ventral <- 
    combined.obj$NESTscore_sample_midbrain <- 
    combined.obj$NESTscore_sample_striatum <- 
    combined.obj$invivo_mixed_dorsal <- 
    combined.obj$invivo_mixed_ventral <- 
    combined.obj$invivo_mixed_midbrain <- 
    combined.obj$invivo_mixed_striatum <- "n/a"
  
  combined.obj$NESTscore_sample_dorsal[rownames(meta_dorsal)] <- meta_dorsal$NESTscore_sample
  combined.obj$NESTscore_sample_ventral[rownames(meta_ventral)] <- meta_ventral$NESTscore_sample
  combined.obj$NESTscore_sample_midbrain[rownames(meta_midbrain)] <- meta_midbrain$NESTscore_sample
  combined.obj$NESTscore_sample_striatum[rownames(meta_striatum)] <- meta_striatum$NESTscore_sample
  
  combined.obj$invivo_mixed_dorsal[rownames(meta_dorsal)] <- meta_dorsal$invivo_mixed
  combined.obj$invivo_mixed_ventral[rownames(meta_ventral)] <- meta_ventral$invivo_mixed
  combined.obj$invivo_mixed_midbrain[rownames(meta_midbrain)] <- meta_midbrain$invivo_mixed
  combined.obj$invivo_mixed_striatum[rownames(meta_striatum)] <- meta_striatum$invivo_mixed

  umap <- data.frame("sample" = combined.obj$sample, 
                     "umap_1" = combined.obj@reductions$umap@cell.embeddings[,1],
                     "umap_2" = combined.obj@reductions$umap@cell.embeddings[,2],
                     "invivo_mixed" = as.character(combined.obj$invivo_mixed),
                     "annotation" = as.character(combined.obj$annotation),
                     "Region" = combined.obj$Region,
                     "protocol_driven" = combined.obj$protocol_driven,
                     "invivo_mixed_dorsal" = combined.obj$invivo_mixed_dorsal,
                     "invivo_mixed_ventral" = combined.obj$invivo_mixed_ventral,
                     "invivo_mixed_midbrain" = combined.obj$invivo_mixed_midbrain,
                     "invivo_mixed_striatum" = combined.obj$invivo_mixed_striatum)
  
  umap$invivo_mixed_dorsal <- c("not well-mixed","well-mixed")[as.numeric(umap$invivo_mixed_dorsal)]
  umap$invivo_mixed_dorsal[umap$invivo_mixed_dorsal == "not well-mixed" & umap$sample == "invivo"] <- "in vivo specific"
  umap$invivo_mixed_dorsal[umap$invivo_mixed_dorsal == "not well-mixed" & umap$sample == "invitro"] <- "in vitro specific"
  
  umap$invivo_mixed_ventral <- c("not well-mixed","well-mixed")[as.numeric(umap$invivo_mixed_ventral)]
  umap$invivo_mixed_ventral[umap$invivo_mixed_ventral == "not well-mixed" & umap$sample == "invivo"] <- "in vivo specific"
  umap$invivo_mixed_ventral[umap$invivo_mixed_ventral == "not well-mixed" & umap$sample == "invitro"] <- "in vitro specific"
  
  umap$invivo_mixed_midbrain <- c("not well-mixed","well-mixed")[as.numeric(umap$invivo_mixed_midbrain)]
  umap$invivo_mixed_midbrain[umap$invivo_mixed_midbrain == "not well-mixed" & umap$sample == "invivo"] <- "in vivo specific"
  umap$invivo_mixed_midbrain[umap$invivo_mixed_midbrain == "not well-mixed" & umap$sample == "invitro"] <- "in vitro specific"
  
  umap$invivo_mixed_striatum <- c("not well-mixed","well-mixed")[as.numeric(umap$invivo_mixed_striatum)]
  umap$invivo_mixed_striatum[umap$invivo_mixed_striatum == "not well-mixed" & umap$sample == "invivo"] <- "in vivo specific"
  umap$invivo_mixed_striatum[umap$invivo_mixed_striatum == "not well-mixed" & umap$sample == "invitro"] <- "in vitro specific"
  
  umap$invivo_mixed[umap$invivo_mixed == "well-mixed"] <- "well-mixed"
  umap$invivo_mixed[umap$invivo_mixed == "not well-mixed" & umap$sample == "invivo"] <- "in vivo specific"
  umap$invivo_mixed[umap$invivo_mixed == "not well-mixed" & umap$sample == "invitro"] <- "in vitro specific"
  
  nr_invitro <- c("telencephalon",
                  "midbrain",
                  "diencephalon",
                  "cerebellum")
  
  umap$Region <- factor(umap$Region, levels = nr_invitro)
  
  saveRDS(umap, "umap_df_integrated_braun_pw14.RDS")
  #saveRDS(combined.obj, "integrated_braun_pw14.RDS")
}
