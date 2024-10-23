rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ragg)
library(patchwork) 
library(pheatmap)
library(NESTscore)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./scRNAseq_analysis/data/"
dir.create("figures/Fig2", showWarnings = FALSE)

combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))
combined.obj$NESTscore_cell.line <- as.numeric(combined.obj$NESTscore_cell.line)

umap_df <- cbind(combined.obj@meta.data, combined.obj@reductions$umap@cell.embeddings)

### Panel 2B ###
agg_png(paste0("figures/Fig2/Fig2B_UMAP_NEST_Score.png"),
    width = 135, height = 90, units = "mm", res = 500, scaling = 0.8)
df_obj <- umap_df[which(combined.obj$protocol == "dorsal"), ]
p1 <- ggplot() +
  geom_point(data = umap_df, aes(x = umap_1, y = umap_2), col = "snow2", size = 0.3, stroke = 0) + 
  geom_point(data = df_obj[order(df_obj$NESTscore_cell.line),], aes(x = umap_1, y = umap_2, color = NESTscore_cell.line), size = 0.3, stroke = 0) + 
  scale_color_viridis_c(name = "NEST-Score") + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) + 
  ggtitle("dorsal protocol")

df_obj <- umap_df[which(combined.obj$protocol == "ventral"), ]
p2 <- ggplot() +
  geom_point(data = umap_df, aes(x = umap_1, y = umap_2), col = "snow2", size = 0.3, stroke = 0) + 
  geom_point(data = df_obj[order(df_obj$NESTscore_cell.line),], aes(x = umap_1, y = umap_2, color = NESTscore_cell.line), size = 0.3, stroke = 0) + 
  scale_color_viridis_c(name = "NEST-Score") + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) + 
  ggtitle("ventral protocol")

df_obj <- umap_df[which(combined.obj$protocol == "midbrain"), ]
p3 <- ggplot() +
  geom_point(data = umap_df, aes(x = umap_1, y = umap_2), col = "snow2", size = 0.3, stroke = 0) + 
  geom_point(data = df_obj[order(df_obj$NESTscore_cell.line),], aes(x = umap_1, y = umap_2, color = NESTscore_cell.line), size = 0.3, stroke = 0) + 
  scale_color_viridis_c(name = "NEST-Score") + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) + 
  ggtitle("midbrain protocol")

df_obj <- umap_df[which(combined.obj$protocol == "striatum"), ]
p4 <- ggplot() +
  geom_point(data = umap_df, aes(x = umap_1, y = umap_2), col = "snow2", size = 0.3, stroke = 0) + 
  geom_point(data = df_obj[order(df_obj$NESTscore_cell.line),], aes(x = umap_1, y = umap_2, color = NESTscore_cell.line), size = 0.3, stroke = 0) + 
  scale_color_viridis_c(name = "NEST-Score") + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black"))  +
  ggtitle("striatum protocol")
print(p1+p2+p3+p4+plot_layout(ncol = 2))
dev.off()

### Panel 2C ###
df_obj <- umap_df[which(combined.obj$protocol == "dorsal"), ]
set.seed(123)

png(paste0("figures/Fig2/Fig2C_UMAP_dorsal_cell.lines.png"),
    width = 45, height = 35, units = "mm", res = 500)
p <- ggplot() +
  geom_point(data = umap_df, aes(x = umap_1, y = umap_2), col = "snow2", size = 0.3, stroke = 0) +
  geom_point(data = df_obj[sample(rownames(df_obj)),], aes(x = umap_1, y = umap_2, color = cell.line), size = 0.3, stroke = 0) +
  scale_color_manual(values = cell.line_colors, name = "cell line") + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=2))) + NoLegend()
print(p)
dev.off()

p <- ggplot() +
  geom_point(data = umap_df, aes(x = umap_1, y = umap_2), col = "snow2", size = 0.3, stroke = 0) +
  geom_point(data = df_obj[sample(rownames(df_obj)),], aes(x = umap_1, y = umap_2, color = cell.line), size = 0.3, stroke = 0) +
  scale_color_manual(values = cell.line_colors, name = "cell line") + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=2))) 
p

pdf("figures/Fig2/Fig2C_UMAP_dorsal_cell.lines_legend.pdf", width = 4, height = 4)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p)))
dev.off()

# cells_dorsal_1 <- CellSelector(DimPlot(combined.obj))
# cells_dorsal_2 <- CellSelector(DimPlot(combined.obj))
# cells_dorsal <- unique(c(cells_dorsal_1, cells_dorsal_2))
# saveRDS(cells_dorsal, paste0(DataDir,"cells_dorsal.RDS"))

cells_dorsal <- readRDS(paste0(DataDir,"cells_dorsal.RDS"))
df_obj_sub <- df_obj[grepl("dorsal",df_obj$annotation_coarse) & 
                       rownames(df_obj) %in% cells_dorsal,]

set.seed(123)
p1 <- ggplot() +
  geom_point(data = df_obj_sub[order(df_obj_sub$NESTscore_cell.line),], aes(x = umap_1, y = umap_2, color = NESTscore_cell.line), 
             size = 0.5, stroke = 0) +
  scale_color_viridis_c(name = "NEST-Score", 
                        limits = range(combined.obj$NESTscore_cell.line[which(combined.obj$protocol == "dorsal")])) + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) 
df_obj_sub <- df_obj_sub[sample(rownames(df_obj_sub)),]
p2 <- ggplot() +
  geom_point(data = df_obj_sub, aes(x = umap_1, y = umap_2, color = cell.line), size = 0.5, stroke = 0) +
  scale_color_manual(values = cell.line_colors, name = "cell line") +
  NoLegend() + theme(legend.position = "none") + theme_void()   + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=2))) #+NoLegend()
try(p1+p2 + ggsave("figures/Fig2/Fig2C_NEST_Score_dorsal_subset_neurons.pdf",
                   width = 100, height = 35, units = "mm"), silent = T)

# cells_stromal <- CellSelector(DimPlot(combined.obj))
# saveRDS(cells_stromal, paste0(DataDir,"cells_stromal.RDS"))

cells_stromal  <- readRDS(paste0(DataDir,"cells_stromal.RDS"))
df_obj_sub <- df_obj[grepl("stromal",df_obj$annotation_coarse) & 
                       rownames(df_obj) %in% cells_stromal,]

set.seed(123)
p1 <- ggplot() +
  geom_point(data = df_obj_sub[order(df_obj_sub$NESTscore_cell.line),], aes(x = umap_1, y = umap_2, color = NESTscore_cell.line), 
             size = 0.5, stroke = 0) +
  scale_color_viridis_c(name = "NEST-Score", 
                        limits = range(combined.obj$NESTscore_cell.line[which(combined.obj$protocol == "dorsal")])) +
  NoLegend() + theme(legend.position = "none")  + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) 

df_obj_sub <- df_obj_sub[sample(rownames(df_obj_sub)),]
p2 <- ggplot() +
  geom_point(data = df_obj_sub, aes(x = umap_1, y = umap_2, color = cell.line), size = 0.5, stroke = 0) +
  scale_color_manual(values = cell.line_colors, name = "cell line") +
  NoLegend() + theme(legend.position = "none") + theme_void()   + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black"))  +
  guides(colour = guide_legend(override.aes = list(size=2)))
try(p1+p2 + ggsave("figures/Fig2/Fig2C_NEST_Score_dorsal_subset_stromal.pdf",
                   width = 100, height = 35, units = "mm"), silent = T)

dev.off()

### Panel 2D ###
pcadim <- 40
k0 <- 100

Idents(combined.obj) <- combined.obj$protocol
ls.protocol <- SplitObject(combined.obj)

for(obj in ls.protocol) {
  
  prot <- as.character(unique(obj$protocol))
  
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:pcadim)
  
  NESTres <- NESTscore(obj, "cell.line", 
                       k_nn = k0, 
                       ndims = pcadim, 
                       return_pairwise_eval = T, 
                       show_heatmap = F, 
                       heatmap_col = obj@misc$protocol_colors[[prot]])
  
  obj <- NESTres$seuratobj
  
  if(prot == "ventral") {
    NESTres$pairwise_matrix <- rbind(cbind(NESTres$pairwise_matrix,"H9" = NA),"H9" = NA)
  }
  
  pdf(paste0("figures/Fig2/Fig2D_Heatmap_cell.lines_relfreq_", prot,".pdf"),
      width = 2.5, height = 2)
  pheatmap::pheatmap(NESTres$pairwise_matrix[levels(obj$cell.line),levels(obj$cell.line)],
                     cluster_rows = F,
                     cluster_cols = F,
                     color = NESTres$proposed_params_pheatmap$color,
                     breaks = NESTres$proposed_params_pheatmap$breaks,
                     legend_breaks = round(seq(0,1,1/length(unique(obj$cell.line))),3),
                     na_col = "transparent",
                     main = paste0(prot," protocol"))
  dev.off()
}
rm(obj,NESTres)
