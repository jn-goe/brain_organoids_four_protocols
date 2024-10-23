rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(dplyr)
library(ggplot2)
library(voxhunt)
load_aba_data('./voxhunt/voxhunt_data')

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./scRNAseq_analysis/data/"
dir.create("figures/Fig1", showWarnings = FALSE)

combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))

round(table(combined.obj$cell.line[combined.obj$annotation_coarse == "stromal" &
                                     combined.obj$protocol == "dorsal"])/
        sum(table(combined.obj$cell.line[combined.obj$annotation_coarse == "stromal" &
                                           combined.obj$protocol == "dorsal"]))*100,2)
# UOFV_1 XUJA_2     H1     H9 ROZH_5    177    178    176 
# 79.31   3.44   0.00   3.61   0.17   0.00   0.22  13.25 

round(table(combined.obj$cell.line[combined.obj$annotation_coarse == "muscle" &
                                     combined.obj$protocol == "dorsal"])/
        sum(table(combined.obj$cell.line[combined.obj$annotation_coarse == "muscle" &
                                           combined.obj$protocol == "dorsal"]))*100,2)
# UOFV_1 XUJA_2     H1     H9 ROZH_5    177    178    176 
# 0.00   4.27   0.02   0.20   0.02   0.00   9.06  86.44 

umap_df <- cbind(combined.obj@meta.data, combined.obj@reductions$umap@cell.embeddings)

### Panel 1B ###
png("figures/Fig1/Fig1B_UMAP_celltypes.png", 
    width = 85, height = 85, units = "mm", res = 500)
p1 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, col = annotation_fine_num)) + 
  geom_point(size = 0.3, stroke = 0, alpha = 1) + scale_color_manual(values = annotation_fine_num_colors) + 
  theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  guides(col = guide_legend(title="", override.aes = list(size = 4))) + NoLegend()
print(p1)
dev.off()

p2 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, col = annotation_fine_num)) + 
  geom_point(size = 0.1, alpha = 1) + scale_color_manual(values = annotation_fine_num_colors) + 
  theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  guides(col = guide_legend(title="",
                             override.aes = list(size = 4))) 

pdf("figures/Fig1/Fig1B_UMAP_celltypes_legend.pdf", width = 5, height = 5)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p2)))
dev.off()

### Panel 1C ###
png("figures/Fig1/Fig1C_UMAP_by_protocol.png", 
    width = 85, height = 85, units = "mm", res = 500)
p1 <- ggplot(umap_df, vars = c("umap_1", "umap_2", "protocol"), 
             aes(x = umap_1, y = umap_2)) +
  ggtitle("dorsal protocol") +
  geom_point(size = 0.2, stroke = 0, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$protocol == "dorsal",], stroke = 0, size = 0.2, alpha=1, col = "seagreen") + 
  theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 
p2 <- ggplot(umap_df, vars = c("umap_1", "umap_2", "protocol"), 
             aes(x = umap_1, y = umap_2)) +
  theme_void() + ggtitle("ventral protocol") +
  geom_point(size = 0.2, stroke = 0, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$protocol == "ventral",], stroke = 0, size = 0.2, alpha=1, col = "firebrick") + 
  theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black"))  
p3 <- ggplot(umap_df, vars = c("umap_1", "umap_2", "protocol"), 
             aes(x = umap_1, y = umap_2)) +
  theme_void() + ggtitle("midbrain protocol") +
  geom_point(size = 0.2, stroke = 0, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$protocol == "midbrain",], stroke = 0, size = 0.2, alpha=1, col = "royalblue4") + 
  theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 
p4 <- ggplot(umap_df, vars = c("umap_1", "umap_2", "protocol"), 
             aes(x = umap_1, y = umap_2)) +
  theme_void() + ggtitle("striatum protocol") +
  geom_point(size = 0.2, stroke = 0, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$protocol == "striatum",], stroke = 0, size = 0.2, alpha=1, col = "gold2") + 
  theme_void() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 
print(cowplot::plot_grid(p1,p2,p3,p4, nrow = 2, ncol = 2))
dev.off()

### Panel 1D ###
p <- ggplot(umap_df, aes(fill = protocol, 
                         x = annotation_coarse)) + 
  geom_bar(position = 'fill') + ylab("ratio of cells") + 
  guides(fill=guide_legend(title="protocol")) + 
  scale_fill_manual(values = protocol_colors) +
  xlab(element_blank()) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black"))
p

try(p + ggsave("figures/Fig1/Fig1D_barplot_celltypes.pdf", 
               width = 110, height = 65, units = "mm"), silent = T) 


### Panel 1E ### 
# voxhunt Analysis #

df <- structure_markers('E13', annotation_level = "custom_2")
regional_markers <- df%>%
  group_by(group) %>%
  top_n(15, auc) %>% 
  {unique(.$gene)}

vox_map <- voxel_map(combined.obj,
                     genes_use = regional_markers, 
                     group_name = "protocol", 
                     stage = "E13", 
                     method = "spearman")

pdf("figures/Fig1/Fig1E_voxhunt_by_protocol.pdf", width = 6, height = 1, pointsize = 8)
(plot_map(vox_map, nrow = 1)[[1]] + theme(text = element_text(size = 8, colour = "black"), 
                                          title = element_text(size = 8, colour = "black"),
                                          legend.text = element_text(size = 8, colour = "black"),
                                          strip.text = element_text(size = 8, colour = "black"))) +
  (plot_map(vox_map, nrow = 1)[[2]] + theme(text = element_text(size = 8, colour = "black"), 
                                            title = element_text(size = 8, colour = "black"),
                                            legend.text = element_text(size = 8, colour = "black"),
                                            strip.text = element_text(size = 8, colour = "black"))) + 
  (plot_map(vox_map, nrow = 1)[[3]] + theme(text = element_text(size = 8, colour = "black"), 
                                            title = element_text(size = 8, colour = "black"),
                                            legend.text = element_text(size = 8, colour = "black"),
                                            strip.text = element_text(size = 8, colour = "black"))) + 
  (plot_map(vox_map, nrow = 1)[[4]] + theme(text = element_text(size = 8, colour = "black"), 
                                            title = element_text(size = 8, colour = "black"),
                                            legend.text = element_text(size = 8, colour = "black"),
                                            strip.text = element_text(size = 8, colour = "black"))) +  
  patchwork::plot_layout(ncol = 4)
dev.off()

pdf("figures/Fig1/Fig1E_voxhunt_reference.pdf", width = 2, height = 1)
voxhunt::plot_annotation('E13', show_legend = T, annotation_level = "custom_1") + theme(text = element_text(size = 8, colour = "black"), 
                                                                                        title = element_text(size = 8, colour = "black"),
                                                                                        legend.text = element_text(size = 8, colour = "black"),
                                                                                        strip.text = element_text(size = 8, colour = "black")) 
dev.off()
