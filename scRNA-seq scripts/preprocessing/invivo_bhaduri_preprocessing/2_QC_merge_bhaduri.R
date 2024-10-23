rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(ggplot2)
source("./helper_functions/mixedness_functions.R")

# Setup ----------------------------
OutDir <- OutDirOrig <- "./scRNAseq_analysis/bhaduri/"
setwd(OutDir)

metadata.bhaduri <- readxl::read_xlsx("41586_2021_3910_MOESM3_ESM.xlsx", sheet = 1)
metadata.bhaduri$cell.name <- sub("gw", "GW", metadata.bhaduri$cell.name)
metadata.bhaduri$cell.name <- paste0(metadata.bhaduri$cell.name,"-1")
metadata.bhaduri <- as.data.frame(metadata.bhaduri)
rownames(metadata.bhaduri) <- metadata.bhaduri$cell.name

set.seed(432)
cells_14 <- metadata.bhaduri[metadata.bhaduri$age == "14" &
                               !(metadata.bhaduri$`cell type` == "Outlier"),]$cell.name
cells_17 <- metadata.bhaduri[metadata.bhaduri$age == "17"
                             & !(metadata.bhaduri$`cell type` == "Outlier"),]$cell.name
cells_18 <- sample_stratified(metadata.bhaduri[metadata.bhaduri$age == "18"  &
                                                 !(metadata.bhaduri$`cell type` == "Outlier"),],
                  group = "cell type",
                  size = 70000 - length(cells_14) - length(cells_17),
                  df_bin = T)
cells_18 <- (metadata.bhaduri[metadata.bhaduri$age == "18",][cells_18,]$cell.name)

OutDir <- OutDirOrig <- "./scRNAseq_analysis/bhaduri/data/"
setwd(OutDir)

objects <- list.files(OutDir)
objects <- objects[grepl(pattern = ".Rds", x = objects)]
objects <- objects[!grepl(pattern = "afterQC", x = objects)]

ls.Seurat <- list()
for(path in rev(objects)) {

  obj <- readRDS(path)
  name <- sub(".Rds","",path)
  print(name)
  obj <- RenameCells(obj, add.cell.id = name)

  if(unique(obj$area) == "occipital") {
    obj$area <- rep("V1", length(obj$area))
    colnames(obj) <- sub("occipital", "V1", colnames(obj))
  }
  if(unique(obj$area) == "somato" & unique(obj$time) == c("GW14")) {
    obj$area <- rep("somatosensory", length(obj$area))
    colnames(obj) <- sub("somato", "somatosensory", colnames(obj))
  }
  if(unique(obj$area) == "ParVZ") {
    obj$area <- rep("ParietalVZ", length(obj$area))
    obj$replicate <- rep("2", length(obj$replicate))
    colnames(obj) <- sub("ParVZ", "ParietalVZ", colnames(obj))
  }
  if(unique(obj$area) == "TempVZ") {
    obj$area <- rep("TemporalVZ", length(obj$area))
    obj$replicate <- rep("2", length(obj$replicate))
    colnames(obj) <- sub("TempVZ", "TemporalVZ", colnames(obj))
  }

  if(unique(obj$time) == "GW14") {
    obj <- subset(obj, cells = cells_14)
  }
  if(unique(obj$time) == "GW17") {
    obj <- subset(obj, cells = cells_17)
  }
  if(unique(obj$time) == "GW18") {
    obj <- subset(obj, cells = cells_18)
  }

  obj$cell_type <- metadata.bhaduri[colnames(obj),]$`cell type`

  obj <- NormalizeData(obj)

  obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^MT\\.|^MT-")
  obj[["percent.mito"]] <- obj[["percent.mito"]]/100

  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPL|^RS")
  obj[["percent.ribo"]] <- obj[["percent.ribo"]]/100

  RPS.genes <- !grepl("^RPL|^RS|^RPS|^MRP", rownames(obj))
  MT.genes <- !grepl("^MT\\.|^MT-", rownames(obj))

  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  VariableFeatures(obj) <- intersect(VariableFeatures(obj),
                                              rownames(obj)[RPS.genes &
                                                                       MT.genes])
  obj <- CellCycleScoring(obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"))

  pca_comp <- min(50, dim(obj)-2)
  obj <- RunPCA(obj, npcs = pca_comp)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:pca_comp, n.neighbors = pca_comp, n.components = 2L)

  ls.Seurat <- append(ls.Seurat, obj)
  saveRDS(obj, paste0(name,"_afterQC.Rds"))
}

saveRDS(ls.Seurat, "bhaduri_GW141718_sampled_object_list.RDS")

combined.obj <- merge(ls.Seurat[[1]],ls.Seurat[2:length(ls.Seurat)])
combined.obj <- JoinLayers(combined.obj)
combined.obj <- AddMetaData(combined.obj, metadata = metadata.bhaduri[colnames(combined.obj),])
ls.Seurat <- SplitObject(combined.obj, split.by = "individual")

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
combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:40)

combined.obj$annotation <- paste0(combined.obj$structure,"_", combined.obj$cell_type)
combined.obj$annotation[combined.obj$cell_type %in% "Endo"] <- "Endo"
combined.obj$annotation[combined.obj$cell_type %in% "Astrocyte"] <- "Astrocyte"
combined.obj$annotation[combined.obj$cell_type %in% "Microglia"] <- "Microglia"
combined.obj$annotation[combined.obj$cell_type %in% "Oligo"] <- "Oligo"
combined.obj$annotation[combined.obj$cell_type %in% "Outlier"] <- "Outlier"
combined.obj$annotation[combined.obj$cell_type %in% "Vascular"] <- "Vascular"
combined.obj$annotation <- tolower(combined.obj$annotation)

combined.obj$annotation <- factor(combined.obj$annotation, levels = c("neocortex_dividing",
                                                                     "neocortex_rg",
                                                                     "neocortex_ipc",
                                                                     "neocortex_neuron",
                                                                     "neocortex_interneuron",

                                                                     "allocortex_dividing",
                                                                     "allocortex_rg",
                                                                     "allocortex_ipc",
                                                                     "allocortex_neuron",
                                                                     "allocortex_interneuron",

                                                                     "claustrum_dividing",
                                                                     "claustrum_rg",
                                                                     "claustrum_interneuron",
                                                                     "claustrum_neuron",

                                                                     "ge_dividing",
                                                                     "ge_rg",
                                                                     "ge_ipc",
                                                                     "ge_neuron",
                                                                     "ge_interneuron",

                                                                     "hypothalamus_dividing",
                                                                     "hypothalamus_rg",
                                                                     "hypothalamus_ipc",
                                                                     "hypothalamus_neuron",

                                                                     "thalamus_dividing",
                                                                     "thalamus_rg",
                                                                     "thalamus_interneuron",
                                                                     "thalamus_neuron",

                                                                     "striatum_dividing",
                                                                     "striatum_rg",
                                                                     "striatum_interneuron",
                                                                     "striatum_neuron",

                                                                     "oligo",
                                                                     "astrocyte",
                                                                     "microglia",

                                                                     "endo",
                                                                     "outlier",
                                                                     "vascular"))

combined.obj@misc$annotation_colors <- annotation_bhaduri_colors

DefaultAssay(combined.obj) <- "RNA"
combined.obj <- JoinLayers(combined.obj)

saveRDS(combined.obj, "bhaduri_GW141718_sampled.RDS")
