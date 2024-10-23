rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(dplyr)
library(ggplot2)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./scRNAseq_analysis/data/"
dir.create("figures/SFig3", showWarnings = FALSE)

invitro <- readRDS(paste0(DataDir, "combined.obj.RDS"))

### SuplFig3A ###
#---- braun in vivo ----#

umap <- readRDS(paste0(DataDir,"umap_df_integrated_braun_pw14.RDS"))
umap_invivo <- cbind(umap[which(umap$sample == "invivo"),])

umap_invitro <- invitro@meta.data
umap_invitro$invivo_mixed <- umap[rownames(umap_invitro),"invivo_mixed"]
umap_invitro$umap_1 <- umap[rownames(umap_invitro),"umap_1"]
umap_invitro$umap_2 <- umap[rownames(umap_invitro),"umap_2"]

#---dorsal---#
round(table(umap_invivo$invivo_mixed_dorsal)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 58.6             41.4 

p1 <- ggplot(umap_invivo, aes(fill=invivo_mixed_dorsal, x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
p1b <- ggplot(umap_invivo, aes(fill=invivo_mixed_dorsal, x="                  all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +     
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff")) + ylab("")

try(cowplot::plot_grid(p1,p1b, rel_widths = c(1,.6)) + 
      ggsave("figures/SFig3/SuplFig3A_barplot_mixed_dorsal.pdf", 
             width = 2, height = 1.5), silent = T)

#---ventral---#
round(table(umap_invivo$invivo_mixed_ventral)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 49.25            50.75 

p1 <- ggplot(umap_invivo, aes(fill=invivo_mixed_ventral, x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
p1b <- ggplot(umap_invivo, aes(fill=invivo_mixed_ventral, x="                  all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff")) + ylab("")
try(cowplot::plot_grid(p1,p1b, rel_widths = c(1,.6)) + 
      ggsave("figures/SFig3/SuplFig3A_barplot_mixed_ventral.pdf", 
             width = 2, height = 1.5), silent = T)

#---midbrain---#
round(table(umap_invivo$invivo_mixed_midbrain)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 72.72            27.28 

p1 <- ggplot(umap_invivo, 
             aes(fill=invivo_mixed_midbrain, x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
p1b <- ggplot(umap_invivo, 
              aes(fill=invivo_mixed_midbrain, x="                  all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff")) + ylab("")

try(cowplot::plot_grid(p1,p1b, rel_widths = c(1,.6)) + 
      ggsave("figures/SFig3/SuplFig3A_barplot_mixed_midbrain.pdf", 
             width = 2, height = 1.5), silent = T)

#---striatum---#
round(table(umap_invivo$invivo_mixed_striatum)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 38.85            61.15 

p1 <- ggplot(umap_invivo, 
             aes(fill=invivo_mixed_striatum, x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
p1b <- ggplot(umap_invivo, 
              aes(fill=invivo_mixed_striatum, x="                  all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff","in vivo specific" = "#d6610aff")) + ylab("")

try(cowplot::plot_grid(p1,p1b, rel_widths = c(1,.6)) + 
      ggsave("figures/SFig3/SuplFig3A_barplot_mixed_striatum.pdf", 
             width = 2, height = 1.5), silent = T)

### SuplFig3B ###

p0 <- ggplot(umap_invitro,  aes(fill=invivo_mixed, x=annotation_coarse)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vitro specific" = "#1b9e76ff"))
try(p0 + ggsave("figures/SFig3/SuplFig3B_Barplot_mixed_celltypes.pdf", width = 3, height = 2), silent = T)

### SuplFig3C ###
#---- bhaduri in vivo ----#

umap <- readRDS(paste0(DataDir,"umap_df_integrated_bhaduri_GW141718_sampled.RDS"))
umap_invivo <- cbind(umap[which(umap$sample == "invivo"),])

umap_invitro <- invitro@meta.data
umap_invitro$invivo_mixed <- umap[rownames(umap_invitro),"invivo_mixed"]
umap_invitro$umap_1 <- umap[rownames(umap_invitro),"umap_1"]
umap_invitro$umap_2 <- umap[rownames(umap_invitro),"umap_2"]

pa <- ggplot(umap_invivo, aes(x=Region)) +
  geom_bar(fill="#d6610aff") + ylab("number of cells") +
  guides(fill=guide_legend(title="protocol")) +
  theme_minimal() + xlab(element_blank()) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black"))
try(pa + ggsave("figures/SFig3/SuplFig3C_cellnumbers_invivo.pdf", 
                width = 2, height = 1.8), silent = T)

### SuplFig3D ###

p0 <- ggplot(umap_invitro, aes(fill=invivo_mixed, x=annotation_coarse)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title=""))  + theme_minimal() +   
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vitro specific" = "#1b9e76ff"))
p0

try(p0 + ggsave("figures/SFig3/SuplFig3D_Barplot_mixed_celltypes.pdf", width = 3, height = 2), silent = T)

### Fig 3E ###
set.seed(123)
p1 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2),
             size = 0.5, stroke = 0,  alpha=1, col = "snow2") +
  geom_point(data = umap_invivo[sample(dim(umap_invivo)[1]),], 
             aes(x = umap_1, y = umap_2, color = annotation),
             size = 0.5, stroke = 0, alpha=1) +
  NoLegend() + 
  scale_color_manual(values = annotation_bhaduri_colors) +
  guides(colour = guide_legend(override.aes = list(size=4), title="")) + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 
try(p1 + ggsave("figures/SFig3/SuplFig3E_UMAP_invivo_celltypes.png", height = 3.5, width = 3.5), silent =T )

set.seed(123)
p1 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2),
             size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_invivo[sample(dim(umap_invivo)[1]),], aes(x = umap_1, y = umap_2, color = annotation),
             size = 0.1, alpha=1) +
  #NoLegend() + 
  scale_color_manual(values = annotation_bhaduri_colors) +
  guides(colour = guide_legend(override.aes = list(size=4), title="")) + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

pdf("figures/SFig3/SuplFig3E_UMAP_invivo_celltypes_legend.pdf", 
    width = 4, height = 4, pointsize = 8)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p1)))
dev.off()

set.seed(123)
p2 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2),
             size = 0.5, alpha=1, stroke = 0, col = "snow2") +
  geom_point(data = umap_invitro[sample(dim(umap_invitro)[1]),], aes(x = umap_1, y = umap_2, color = annotation_fine_num),
             size = 0.5, alpha=1, stroke = 0) + NoLegend() + 
  scale_color_manual(values = invitro@misc$annotation_fine_num_colors) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 
try(p2 + NoLegend() + ggsave("figures/SFig3/SuplFig3E_UMAP_invitro_celltypes.png", height = 3.5, width = 3.5), silent =T )

set.seed(123)
p2 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2),
             size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_invitro[sample(dim(umap_invitro)[1]),], aes(x = umap_1, y = umap_2, color = annotation_fine_num),
             size = 0.1, alpha=.1) + 
  scale_color_manual(values = invitro@misc$annotation_fine_num_colors) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

pdf("figures/SFig3/SuplFig3E_UMAP_invitro_celltypes_legend.pdf", width = 4, height = 4, pointsize = 8)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p2)))
dev.off()

### Fig 3F ###
p0 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2),
             size = 0.5, alpha=1, stroke = 0, col = "snow2") +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2, color = invivo_mixed),
             size = 0.5, alpha=1, stroke = 0) +
  scale_color_manual(values = c("in vivo specific" = "#d6610aff", 
                                "in vitro specific" = "#1b9e76ff",
                                "well-mixed" = "#6964aaff")) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend()
try(p0 + ggsave("figures/SFig3/SuplFig3F_UMAP_Mix_Score_invitro.png", height = 4, width = 4), silent =T )

p0 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2),
             size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2, color = invivo_mixed),
             size = 0.1, alpha=1) +
  scale_color_manual(values = c("in vivo specific" = "#d6610aff", 
                                "in vitro specific" = "#1b9e76ff",
                                "well-mixed" = "#6964aaff")) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

pdf("figures/SFig3/SuplFig3F_UMAP_Mix_Score_invitro_legend.pdf", width = 4, height = 4, pointsize = 8)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p0)))
dev.off()

p2 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2),
             size = 0.5, alpha=1, stroke = 0, col = "snow2") +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2, color = invivo_mixed),
             size = 0.5, alpha=1, stroke = 0) +
  scale_color_manual(values = c("in vivo specific" = "#d6610aff", 
                                "in vitro specific" = "#1b9e76ff",
                                "well-mixed" = "#6964aaff")) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend()
try(p2 + ggsave("figures/SFig3/SuplFig3F_UMAP_Mix_Score_invivo.png", height = 4, width = 4), silent =T )

p2 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2),
             size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2, color = invivo_mixed),
             size = 0.1, alpha=1) +
  scale_color_manual(values = c("in vivo specific" = "#d6610aff", 
                                "in vitro specific" = "#1b9e76ff",
                                "well-mixed" = "#6964aaff")) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

pdf("figures/SFig3/SuplFig3F_UMAP_Mix_Score_invivo_legend.pdf", width = 4, height = 4, pointsize = 8)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p2)))
dev.off()

#### Fig 3G ###
round(table(umap_invivo$invivo_mixed)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 22.96            77.04 

p0 <- ggplot(umap_invivo, aes(fill=invivo_mixed, 
                                                                           x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
p0b <- ggplot(umap_invivo, aes(fill=invivo_mixed, 
                                                                            x="                              all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff")) + ylab("")

try(cowplot::plot_grid(p0,p0b, rel_widths = c(.5,.3), rel_heights = ) + 
      ggsave("figures/SFig3/SuplFig3G_Mix_Score_invivo_cellnumbers.pdf", 
             height = 2.5, width = 2), silent =T )

### SuplFig3H ###
#---dorsal---#
round(table(umap_invivo$invivo_mixed_dorsal)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 56.29            43.71 

p1 <- ggplot(umap_invivo, 
             aes(fill=invivo_mixed_dorsal, x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
p1b <- ggplot(umap_invivo, 
              aes(fill=invivo_mixed_dorsal, x="                              all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +     
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff")) + ylab("")

try(cowplot::plot_grid(p1,p1b, rel_widths = c(1,.6)) + 
      ggsave("figures/SFig3/SuplFig3H_barplot_mixed_dorsal.pdf", 
             width = 2, height = 1.5), silent = T)

#---ventral---#
round(table(umap_invivo$invivo_mixed_ventral)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 50.6             49.4 

p1 <- ggplot(umap_invivo, 
             aes(fill=invivo_mixed_ventral, x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
p1b <- ggplot(umap_invivo, 
              aes(fill=invivo_mixed_ventral, x="                              all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff")) + ylab("")

try(cowplot::plot_grid(p1,p1b, rel_widths = c(1,.6)) + 
      ggsave("figures/SFig3/SuplFig3H_barplot_mixed_ventral.pdf", 
             width = 2, height = 1.5), silent = T)

#---midbrain---#
round(table(umap_invivo$invivo_mixed_midbrain)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 67.75            32.25 

p1 <- ggplot(umap_invivo, 
             aes(fill=invivo_mixed_midbrain, x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
#p1
p1b <- ggplot(umap_invivo, 
              aes(fill=invivo_mixed_midbrain, x="                              all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff")) + ylab("")

try(cowplot::plot_grid(p1,p1b, rel_widths = c(1,.6)) + 
      ggsave("figures/SFig3/SuplFig3H_barplot_mixed_midbrain.pdf", 
             width = 2, height = 1.5), silent = T)

#---striatum---#
round(table(umap_invivo$invivo_mixed_striatum)/dim(umap_invivo)[1]*100,2)
# in vivo specific       well-mixed 
# 29.91            70.09 

p1 <- ggplot(umap_invivo, 
             aes(fill=invivo_mixed_striatum, x=Region)) + 
  geom_bar(position="fill") + ylab("ratio") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff"))
p1b <- ggplot(umap_invivo, 
              aes(fill=invivo_mixed_striatum, x="                              all")) + 
  geom_bar(position="fill") + ylab("") + xlab("") +
  guides(fill=guide_legend(title="")) + theme_minimal() +       
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1.05),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend() + 
  scale_fill_manual(values = c("well-mixed" = "#6964aaff",
                               "in vivo specific" = "#d6610aff")) + ylab("")

try(cowplot::plot_grid(p1,p1b, rel_widths = c(1,.6)) + 
      ggsave("figures/SFig3/SuplFig3H_barplot_mixed_striatum.pdf", 
             width = 2, height = 1.5), silent = T)
