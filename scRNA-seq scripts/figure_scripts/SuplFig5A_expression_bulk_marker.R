rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(patchwork)

source("./helper_functions/mixedness_functions.R")
source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
OutDir <- OutDirOrig <- "./"
setwd(OutDirOrig)
dir.create("figures/SFig5", showWarnings = FALSE)

DataDirsc <- "./scRNAseq_analysis/data/"
DataDirbulk <- "./bulkRNAseq_analysis/data/"

combined.obj <- readRDS(paste0(DataDirsc,"combined.obj.RDS"))
marker <- read.csv(paste0(DataDirbulk,"bulk_RNAseq_marker.csv"))

marker$X <- NULL
marker_dorsal <- marker[marker$protocol == "dorsal", ][1:50,]
marker_ventral <- marker[marker$protocol == "ventral", ][1:50,]
marker_midbrain <- marker[marker$protocol == "midbrain", ][1:50,]
marker_striatum <- marker[marker$protocol == "striatum", ][1:50,]

combined.obj <- AddModuleScore(combined.obj, list(marker_dorsal$Gene), name = "dorsal")
combined.obj <- AddModuleScore(combined.obj, list(marker_ventral$Gene), name = "ventral")
combined.obj <- AddModuleScore(combined.obj, list(marker_midbrain$Gene), name = "midbrain")
combined.obj <- AddModuleScore(combined.obj, list(marker_striatum$Gene), name = "striatum")

combined.obj$bulk_marker_dorsal1 <- combined.obj$dorsal1
combined.obj$bulk_marker_dorsal1[combined.obj$bulk_marker_dorsal1 > quantile(combined.obj$bulk_marker_dorsal1, 0.99)] <- quantile(combined.obj$bulk_marker_dorsal1, 0.99)
combined.obj$bulk_marker_dorsal1[combined.obj$bulk_marker_dorsal1 < quantile(combined.obj$bulk_marker_dorsal1, 0.01)] <- quantile(combined.obj$bulk_marker_dorsal1, 0.01)

combined.obj$bulk_marker_ventral1 <- combined.obj$ventral1
combined.obj$bulk_marker_ventral1[combined.obj$bulk_marker_ventral1 > quantile(combined.obj$bulk_marker_ventral1, 0.99)] <- quantile(combined.obj$bulk_marker_ventral1, 0.99)
combined.obj$bulk_marker_ventral1[combined.obj$bulk_marker_ventral1 < quantile(combined.obj$bulk_marker_ventral1, 0.01)] <- quantile(combined.obj$bulk_marker_ventral1, 0.01)

combined.obj$bulk_marker_midbrain1 <- combined.obj$midbrain1
combined.obj$bulk_marker_midbrain1[combined.obj$bulk_marker_midbrain1 > quantile(combined.obj$bulk_marker_midbrain1, 0.99)] <- quantile(combined.obj$bulk_marker_midbrain1, 0.99)
combined.obj$bulk_marker_midbrain1[combined.obj$bulk_marker_midbrain1 < quantile(combined.obj$bulk_marker_midbrain1, 0.01)] <- quantile(combined.obj$bulk_marker_midbrain1, 0.01)

combined.obj$bulk_marker_striatum1 <- combined.obj$striatum1
combined.obj$bulk_marker_striatum1[combined.obj$bulk_marker_striatum1 > quantile(combined.obj$bulk_marker_striatum1, 0.99)] <- quantile(combined.obj$bulk_marker_striatum1, 0.99)
combined.obj$bulk_marker_striatum1[combined.obj$bulk_marker_striatum1 < quantile(combined.obj$bulk_marker_striatum1, 0.01)] <- quantile(combined.obj$bulk_marker_striatum1, 0.01)

df <- cbind(combined.obj@reductions$umap@cell.embeddings, combined.obj@meta.data)

df_dorsal <- df
df_dorsal <- df_dorsal[order(df_dorsal$bulk_marker_dorsal1),]
p1a <- ggplot(df_dorsal, aes(x =umap_1, y = umap_2, col = bulk_marker_dorsal1)) + 
  geom_point(data = df, aes(x = , y = umap_2), col = "snow2", size = 0.5, stroke = 0) + 
  scale_color_viridis_c(name = "") + 
  geom_point(size = 0.5, stroke = 0) + ggtitle("dorsal protocol marker") + theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

df_ventral <- df
df_ventral <- df_ventral[order(df_ventral$bulk_marker_ventral1),]
p2a <- ggplot(df_ventral, aes(x =umap_1, y = umap_2, col = bulk_marker_ventral1)) + 
  scale_color_viridis_c(name = "") + 
  geom_point(data = df, aes(x = , y = umap_2), col = "snow2", size = 0.5, stroke = 0) + 
  geom_point(size = 0.5, stroke = 0) + ggtitle("ventral protocol marker") + theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

df_midbrain <- df
df_midbrain <- df_midbrain[order(df_midbrain$bulk_marker_midbrain1),]
p3a <- ggplot(df_midbrain, aes(x =umap_1, y = umap_2, col = bulk_marker_midbrain1)) + 
  scale_color_viridis_c(name = "") + 
  geom_point(data = df, aes(x = , y = umap_2), col = "snow2", size = 0.5, stroke = 0) + 
  geom_point(size = 0.5, stroke = 0) + ggtitle("midbrain protocol marker") + theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

df_striatum <- df
df_striatum <- df_striatum[order(df_striatum$bulk_marker_striatum1),]
p4a <- ggplot(df_striatum, aes(x = umap_1, y = umap_2, col = bulk_marker_striatum1)) + 
  scale_color_viridis_c(name = "") + 
  geom_point(data = df, aes(x = , y = umap_2), col = "snow2", size = 0.5, stroke = 0) + 
  geom_point(size = 0.5, stroke = 0) + ggtitle("striatum protocol marker") + theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

try(p1a+p2a+p3a+p4a + plot_layout(nrow = 1) + 
      ggsave("figures/SFig5/SuplFig5A_UMAP_bulk_marker.png", width = 10, height = 2),silent = T)
