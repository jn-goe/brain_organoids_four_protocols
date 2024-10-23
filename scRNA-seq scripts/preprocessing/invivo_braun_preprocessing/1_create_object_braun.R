rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(anndata)

# Setup ----------------------------
setwd("./scRNAseq_analysis/linnarsson/")

###### transform from AnnData to Seurat ###### 
genes <- read.csv("genes.csv")
gene_names <- sub("'","",sub("b'","",genes$Gene)) 

f <- "pw14.h5ad"

adata <- read_h5ad(f, backed = T)
adata <- t(as.matrix(adata$X))
rownames(adata) <- make.unique(gene_names)

adata <- adata[genes$ValidGenes == "True", ]

meta <- read.csv(paste0("cellmetadata",sub(".h5ad",".csv",sub(".","",sub("p","",f)))))
for(c in colnames(meta)) { 
  meta[,c] <- sub("'","",sub("b'","",meta[,c]))
}
rownames(meta) <- meta$X
meta$X <- NULL
colnames(adata) <- rownames(meta)

HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(adata), 
                                             unmapped.as.na = FALSE, 
                                             map = NULL, 
                                             species = "human")
rownames(adata) <- make.unique(HGNC.updated$Suggested.Symbol)

obj <- CreateSeuratObject(adata, min.cell = 10)
obj <- AddMetaData(obj, meta)

RPS.genes <- !grepl("^RPL|^RS|^RPS|^MRP", rownames(obj))
MT.genes <- !grepl("^MT\\.|^MT-", rownames(obj))
MALAT1.genes <- !grepl("MALAT1", rownames(obj))

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
VariableFeatures(obj) <- intersect(VariableFeatures(obj),
                                            rownames(obj)[RPS.genes & MT.genes & MALAT1.genes])
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:40, n.components = 2L)

obj$Region <- tolower(obj$Region)
obj$Region <- factor(obj$Region, levels = c("telencephalon",
                                                              "midbrain",
                                                              "diencephalon",
                                                              "cerebellum"))

obj$Subregion <- tolower(obj$Subregion)
obj$Subregion[which(obj$Subregion == "midbrain dorsal")] <- "dorsal midbrain"
obj$Subregion[which(obj$Subregion == "midbrain ventral")] <- "ventral midbrain"
obj$Subregion <- factor(obj$Subregion, levels = c("cortex",
                                                                    "hippocampus",     
                                                                    "striatum",         
                                                                    "dorsal midbrain",  
                                                                    "ventral midbrain",
                                                                    "thalamus",
                                                                    "hypothalamus",     
                                                                    "cerebellum"))

obj$CellClass <- tolower(obj$CellClass)
obj$CellClass <- factor(obj$CellClass, levels = c("neural crest", 
                                                                    "glioblast",
                                                                    "radial glia",  
                                                                    "neuroblast",  
                                                                    "neuronal ipc", 
                                                                    "neuron",       
                                                                    "oligo",        
                                                                    "fibroblast",   
                                                                    "immune",      
                                                                    "vascular",
                                                                    "erythrocyte"))

obj$annotation <- paste0(obj$Subregion, "_", obj$CellClass)
obj$annotation[obj$CellClass == "neural crest"] <- "neural crest"
obj$annotation[obj$CellClass == "oligo"] <- "oligo"
obj$annotation[obj$CellClass == "fibroblast"] <- "fibroblast"
obj$annotation[obj$CellClass == "immune"] <- "immune"
obj$annotation[obj$CellClass == "vascular"] <- "vascular"
obj$annotation[obj$CellClass == "erythrocyte"] <- "erythrocyte"

obj$annotation <- factor(obj$annotation, levels = c("neural crest",
                                                                  
                                                                  "cortex_glioblast",
                                                                  "cortex_radial glia",
                                                                  "cortex_neuroblast",
                                                                  "cortex_neuronal ipc",
                                                                  "cortex_neuron",
                                                                  
                                                                  "hippocampus_glioblast",
                                                                  "hippocampus_radial glia",
                                                                  "hippocampus_neuroblast",
                                                                  "hippocampus_neuronal ipc",
                                                                  "hippocampus_neuron",
                                                                  
                                                                  "striatum_glioblast",
                                                                  "striatum_radial glia",
                                                                  "striatum_neuroblast",
                                                                  "striatum_neuronal ipc",
                                                                  "striatum_neuron",
                                                                  
                                                                  "dorsal midbrain_glioblast",
                                                                  "dorsal midbrain_radial glia",
                                                                  "dorsal midbrain_neuroblast",
                                                                  "dorsal midbrain_neuronal ipc",
                                                                  "dorsal midbrain_neuron",
                                                                  
                                                                  "ventral midbrain_glioblast",
                                                                  "ventral midbrain_radial glia",
                                                                  "ventral midbrain_neuroblast",
                                                                  "ventral midbrain_neuronal ipc",
                                                                  "ventral midbrain_neuron",
                                                                  
                                                                  "thalamus_glioblast",
                                                                  "thalamus_radial glia",
                                                                  "thalamus_neuroblast",
                                                                  "thalamus_neuronal ipc",
                                                                  "thalamus_neuron",
                                                                  
                                                                  "hypothalamus_glioblast",
                                                                  "hypothalamus_radial glia",
                                                                  "hypothalamus_neuroblast",
                                                                  "hypothalamus_neuronal ipc",
                                                                  "hypothalamus_neuron",
                                                                  
                                                                  "cerebellum_glioblast",
                                                                  "cerebellum_radial glia",
                                                                  "cerebellum_neuroblast",
                                                                  "cerebellum_neuronal ipc",
                                                                  "cerebellum_neuron",
                                                                  
                                                                  "oligo",        
                                                                  "fibroblast",   
                                                                  "immune",      
                                                                  "vascular",
                                                                  "erythrocyte"))

obj@misc$Region_colors <- Region_colors
obj@misc$Subregion_colors <- Subregion_colors
obj@misc$CellClass_colors <- CellClass_colors
obj@misc$annotation_colors <- annotation_braun_colors

saveRDS(obj, "braun.obj.processed.14.RDS")
