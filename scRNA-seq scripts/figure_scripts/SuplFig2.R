rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(dplyr)
library(ggplot2)
library(NESTscore)
library(GO.db)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./scRNAseq_analysis/data/"
dir.create("figures/SFig2", showWarnings = FALSE)

combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))
combined.obj$NESTscore_cell.line <- as.numeric(combined.obj$NESTscore_cell.line)

umap_df <- cbind(combined.obj@meta.data, combined.obj@reductions$umap@cell.embeddings)

### Panel 2A ###
pcadim <- 40
k0 <- 100

Idents(combined.obj) <- combined.obj$protocol
ls.protocol <- SplitObject(combined.obj)
obj <- ls.protocol[["dorsal"]]
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:pcadim)

plot_thresh <- NEST_search_k(obj, group_by = "cell.line", 
                             k_search_vector = seq(2,200), 
                             k_initial = 100, 
                             thresh_initial = combined.obj@misc$protocol_driven[["dorsal"]])

try((plot_thresh+NoLegend()) + ggsave(paste0("figures/SFig2/SuplFig2A_NEST_Score_treshold_dorsal.png"), width = 2.5, height = 2.5), silent = T)

pdf("figures/SFig2/SuplFig2A_NEST_Score_treshold_dorsal_legend.pdf", width = 2, height = 5, pointsize = 8)
print(grid::grid.newpage() +
        grid::grid.draw(cowplot::get_legend(plot_thresh)))
dev.off()

### Panel 2B ###
combined.obj$protocol_driven_plot <- paste0(combined.obj$protocol, " ", combined.obj$protocol_driven)
combined.obj$protocol_driven_plot[grepl("cell line-driven",combined.obj$protocol_driven_plot)] <- "cell line-driven"
combined.obj$protocol_driven_plot <- factor(combined.obj$protocol_driven_plot, levels = c("cell line-driven", 
                                                                                          "dorsal protocol-driven",
                                                                                          "ventral protocol-driven",
                                                                                          "midbrain protocol-driven",
                                                                                          "striatum protocol-driven"))
p <- ggplot() + xlab("") + ylab("ratio") +
  geom_bar(data = combined.obj@meta.data, 
           aes(fill = protocol_driven_plot,
               x = cell.line), position = "fill") + 
  guides(fill=guide_legend(title="")) + 
  scale_fill_manual(values = c("cell line-driven" = "grey70", 
                               "dorsal protocol-driven" = protocol_colors[["dorsal"]],
                               "ventral protocol-driven" = protocol_colors[["ventral"]],
                               "midbrain protocol-driven" = protocol_colors[["midbrain"]],
                               "striatum protocol-driven" = protocol_colors[["striatum"]])) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) +
  facet_wrap(~protocol)

try(p + ggsave("figures/SFig2/SuplFig2B_barplot_protocol_driven_per_cell.line.pdf", 
               width = 100, height = 60, units = "mm"), silent = T)

### Panel 2C ###
thresh_df <- data.frame("protocol" = names(combined.obj@misc$protocol_driven),
                        "thresh" = combined.obj@misc$protocol_driven)
thresh_df$protocol <- factor(thresh_df$protocol, levels = names(protocol_colors[1:4]))

p <- ggplot(data = combined.obj@meta.data, 
            aes(fill = annotation_coarse,
                x = annotation_coarse, 
                y = NESTscore_cell.line)) + 
  xlab("") + ylab("NEST-Score") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5) + 
  guides(fill=guide_legend(title="")) + 
  scale_fill_manual(values = unlist(combined.obj@misc$annotation_coarse_colors)) + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) +
  geom_abline(data = thresh_df,
              aes(intercept = thresh, slope = 0), lty = 2) +
  facet_wrap(~protocol) + NoLegend()

try(p + ggsave("figures/SFig2/SuplFig2C_boxplot_NEST_Score_per_celltype.pdf", 
               width = 190, height = 100, units = "mm"), silent = T)

### Panel 2D ###

# load Amin et al. data 
obj <- readRDS("./GSE233574_OrganoidScreen_processed_SeuratObject.rds")
pcadim <- 50
k0 <- 100

NESTres <- NESTscore(obj, "condition", 
                     k_nn = k0, 
                     ndims = pcadim,
                     return_pairwise_eval = T, 
                     show_heatmap = F,
                     scale_column = F)
obj <- NESTres$seuratobj

umap_df <- cbind(obj@meta.data, obj@reductions$umap@cell.embeddings)

p1 <- ggplot() +
  geom_point(data = umap_df,aes(x = UMAP_1, y = UMAP_2, color = NESTscore_condition), size = 0.3, stroke = 0) + 
  scale_color_viridis_c(name = "NEST-Score") + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) + 
  ggtitle("") + scPhylo::guide_small_legend_cont() 
ggsave(filename = "figures/SFig2/SuplFig2D_Aminetal_UMAP_NEST_Score.png", plot = p1, width = 60, height = 51, units = "mm", dpi = 500)

p2 <- ggplot() +
  geom_point(data = umap_df,aes(x = UMAP_1, y = UMAP_2, color = region), size = 0.3, stroke = 0) + 
  theme_void() + theme(text = element_text(size = 8, colour = "black"), 
                       title = element_text(size = 8, colour = "black"),
                       legend.text = element_text(size = 8, colour = "black"),
                       strip.text = element_text(size = 8, colour = "black")) + 
  ggtitle("") + scPhylo::guide_small_legend_disc() + 
  scale_color_discrete(name = "region")
ggsave(filename = "figures/SFig2/SuplFig2D_Aminetal_UMAP_region_class.png", p2, width = 65, height = 51, units = "mm", dpi = 500)

### Panel 2E ###

png("figures/SFig2/SuplFig2E_Aminetal_Heatmap_morphogens.png",units = "in", res = 900,
    width = 4.5, height = 3.6)
pheatmap::pheatmap(NESTres$pairwise_matrix, 
                   cluster_rows = T,
                   cluster_cols = T,
                   color = NESTres$proposed_params_pheatmap$color,
                   breaks = NESTres$proposed_params_pheatmap$breaks,
                   clustering_distance_cols = "correlation",
                   clustering_distance_rows = "correlation",
                   na_col = "transparent", 
                   treeheight_col = 5, 
                   treeheight_row = 3,
                   fontsize = 6)
dev.off()
