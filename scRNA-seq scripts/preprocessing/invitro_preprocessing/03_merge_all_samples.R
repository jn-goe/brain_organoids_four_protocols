rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(ggplot2)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
OutDirOrig <- "./scRNAseq_analysis/separate_organoids/"
setwd(OutDirOrig)

DataDir <- "./scRNAseq_analysis/data/"

folderlist <- list.files(OutDirOrig)
folderlist <- folderlist[grepl("organoids", folderlist)]

ls.Seurat <- list()
for(f_index in folderlist) {
  obj <- readRDS(paste0(f_index,"/combined.obj.RDS"))
  ls.Seurat <- append(ls.Seurat, obj)
}
combined.obj <- merge(ls.Seurat[[1]], ls.Seurat[2:length(ls.Seurat)])

combined.obj$protocol <- as.character(combined.obj$protocol)
combined.obj$protocol[combined.obj$protocol == "D"] <- "dorsal"
combined.obj$protocol[combined.obj$protocol == "V"] <- "ventral"
combined.obj$protocol[combined.obj$protocol == "M"] <- "midbrain"
combined.obj$protocol[combined.obj$protocol == "S"] <- "striatum"

combined.obj$protocol <- factor(combined.obj$protocol, 
                                levels = names(protocol_colors)[-length(names(protocol_colors))])

combined.obj$cell.line.repl <- as.character(combined.obj$cellline)
combined.obj$cell.line <- as.character(combined.obj$cellline)
combined.obj$cell.line[combined.obj$cell.line %in% c("176_E1","176_E2")] <- "176"
combined.obj$cellline <- NULL

combined.obj$cell.line.repl <- factor(combined.obj$cell.line.repl, 
                                levels = names(cell.line.repl_colors))
combined.obj$cell.line <- factor(combined.obj$cell.line, 
                                 levels = names(cell.line_colors))

rm(ls.Seurat, obj)

combined.obj <- JoinLayers(combined.obj)

new_counts <- combined.obj@assays$RNA@layers$counts
colnames(new_counts) <- colnames(combined.obj)
rownames(new_counts) <- rownames(combined.obj)
combined.obj <- CreateSeuratObject(new_counts, 
                                   min.cell = 10, 
                                   meta.data = combined.obj@meta.data)

combined.obj@misc$cell.line_colors <- cell.line_colors
combined.obj@misc$cell.line.repl_colors <- cell.line.repl_colors
combined.obj@misc$protocol_colors <- protocol_colors
combined.obj@misc$protocol_colors_bright <- protocol_colors_bright

rm(new_counts)

OutDirOrig <- "./scRNAseq_analysis/"
OutDir <- paste0(OutDirOrig, "allinone/")
setwd(OutDir)

RPS.genes <- !grepl("^RPL|^RS|^RPS|^MRP", rownames(combined.obj))
MT.genes <- !grepl("^MT\\.|^MT-", rownames(combined.obj))

combined.obj <- NormalizeData(combined.obj)
combined.obj <- CellCycleScoring(combined.obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
combined.obj <- FindVariableFeatures(combined.obj, nfeatures = 3000)
VariableFeatures(combined.obj) <- intersect(VariableFeatures(combined.obj),
                                            rownames(combined.obj)[RPS.genes &
                                                                     MT.genes])
combined.obj <- ScaleData(combined.obj, vars.to.regress = c("S.Score", "G2M.Score"))
combined.obj <- RunPCA(combined.obj)
combined.obj <- RunUMAP(combined.obj, dims = 1:40, n.components = 2L, seed.use = 5) 
combined.obj <- FindNeighbors(combined.obj, dims = 1:40, reduction = "pca")

saveRDS(combined.obj, paste0(DataDir,"combined.obj.RDS"))
