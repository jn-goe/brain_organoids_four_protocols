rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./scRNAseq_analysis/data/"
dir.create("figures/Fig3", showWarnings = FALSE)

invitro <- readRDS(paste0(DataDir, "combined.obj.RDS"))

#---- braun in vivo ----#
umap <- readRDS(paste0(DataDir,"umap_df_integrated_braun_pw14.RDS"))
umap_invivo <- cbind(umap[which(umap$sample == "invivo"),])

umap_invitro <- invitro@meta.data
umap_invitro$invivo_mixed <- umap[rownames(umap_invitro),"invivo_mixed"]
umap_invitro$umap_1 <- umap[rownames(umap_invitro),"umap_1"]
umap_invitro$umap_2 <- umap[rownames(umap_invitro),"umap_2"]

### Fig 3A ###
pa <- ggplot(umap_invivo, aes(x=Region)) +
  geom_bar(fill="#d6610aff") + ylab("number of cells") +
  guides(fill=guide_legend(title="protocol")) +
  theme_minimal() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  xlab(element_blank()) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
try(pa + ggsave("figures/Fig3/Fig3A_cellnumbers_invivo.pdf", width = 40, height = 35, units = "mm"), silent = T)

pb <- ggplot(invitro@meta.data, aes(x=protocol)) +
  geom_bar(fill="#1b9e76ff") + ylab("number of cells") +
  guides(fill=guide_legend(title="protocol")) +
  theme_minimal() + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  xlab(element_blank()) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
try(pb + ggsave("figures/Fig3/Fig3A_cellnumbers_invitro.pdf", width = 40, height = 35, units = "mm"), silent = T)

### Fig 3B ###
set.seed(123)
p1 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2),
             size = 0.3, alpha=1, stroke = 0, col = "snow2") +
  geom_point(data = umap_invivo[sample(dim(umap_invivo)[1]),], aes(x = umap_1, y = umap_2, color = annotation),
             size = 0.3, alpha=1, stroke = 0) +
  NoLegend() + 
  scale_color_manual(values = annotation_braun_colors) +
  guides(colour = guide_legend(override.aes = list(size=4), title="")) + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 
try(p1 + ggsave("figures/Fig3/Fig3B_UMAP_invivo_celltypes.png", height = 3.5, width = 3.5), silent =T )

set.seed(123)
p1 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2),
             size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_invivo[sample(dim(umap_invivo)[1]),], aes(x = umap_1, y = umap_2, color = annotation),
             size = 0.1, alpha=1) +
  #NoLegend() + 
  scale_color_manual(values = annotation_braun_colors) +
  guides(colour = guide_legend(override.aes = list(size=4), title="")) + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

pdf("figures/Fig3/Fig3B_UMAP_invivo_celltypes_legend.pdf", 
    width = 3.5, height = 3.5, pointsize = 8)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p1)))
dev.off()

set.seed(123)
p2 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2),
             size = 0.3, alpha=1, stroke = 0, col = "snow2") +
  geom_point(data = umap_invitro[sample(dim(umap_invitro)[1]),], aes(x = umap_1, y = umap_2, color = annotation_fine_num),
             size = 0.3, alpha=1, stroke = 0) + NoLegend() + 
  scale_color_manual(values = invitro@misc$annotation_fine_num_colors) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 
try(p2 + NoLegend() + ggsave("figures/Fig3/Fig3B_UMAP_invitro_celltypes.png", height = 3.5, width = 3.5), silent =T )

set.seed(123)
p2 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2),
             size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_invitro[sample(dim(umap_invitro)[1]),], aes(x = umap_1, y = umap_2, color = annotation_fine_num),
             size = 0.1, alpha=1) + 
  scale_color_manual(values = invitro@misc$annotation_fine_num_colors) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) 

pdf("figures/Fig3/Fig3B_UMAP_invitro_celltypes_legend.pdf", width = 3.5, height = 3.5, pointsize = 8)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p2)))
dev.off()

### Fig 3C ###

p0 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2),
             size = 0.3, alpha=1, stroke = 0, col = "snow2") +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2, color = invivo_mixed),
             size = 0.3, alpha=1, stroke = 0) +
  scale_color_manual(values = c("in vivo specific" = "#d6610aff", 
                                "in vitro specific" = "#1b9e76ff",
                                "well-mixed" = "#6964aaff")) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend()
try(p0 + ggsave("figures/Fig3/Fig3C_UMAP_NEST_Score_invitro.png", height = 4, width = 4), silent =T )

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

pdf("figures/Fig3/Fig3C_UMAP_NEST_Score_invitro_legend.pdf", width = 3.5, height = 3.5, pointsize = 8)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p0)))
dev.off()

p2 <- ggplot() +
  theme_void() +
  geom_point(data = umap_invitro, aes(x = umap_1, y = umap_2),
             size = 0.3, alpha=1, stroke = 0, col = "snow2") +
  geom_point(data = umap_invivo, aes(x = umap_1, y = umap_2, color = invivo_mixed),
             size = 0.3, alpha=1, stroke = 0) +
  scale_color_manual(values = c("in vivo specific" = "#d6610aff", 
                                "in vitro specific" = "#1b9e76ff",
                                "well-mixed" = "#6964aaff")) +
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  + 
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + NoLegend()
try(p2 + ggsave("figures/Fig3/Fig3C_UMAP_NEST_Score_invivo.png", height = 4, width = 4), silent =T )

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

pdf("figures/Fig3/Fig3C_UMAP_NEST_Score_invivo_legend.pdf", width = 3.5, height = 3.5, pointsize = 8)
print(grid::grid.newpage() + 
        grid::grid.draw(cowplot::get_legend(p2)))
dev.off()

#### Fig 3D ###
round(table(umap$invivo_mixed[umap$sample == "invivo"])/
        sum(umap$sample == "invivo")*100, 2)
# in vivo specific       well-mixed 
# 29.39            70.61 

p0 <- ggplot(umap_invivo, aes(fill=invivo_mixed, x=Region)) + 
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
p0b <- ggplot(umap_invivo, aes(fill=invivo_mixed, x="                   all")) + 
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

try(cowplot::plot_grid(p0,p0b, rel_widths = c(.5,.3)) + 
      ggsave("figures/Fig3/Fig3D_NEST_Score_invivo_cellnumbers.pdf", 
             height = 2.5, width = 2), silent =T )

### Fig3E ###

umap$combi <- paste0(umap$protocol_driven, " - ", umap$invivo_mixed)
umap$combi[umap$sample == "invivo"] <- NA

invitro$combi <- umap[colnames(invitro), "combi"]

umap_df <- as.data.frame(cbind(invitro@meta.data, invitro@reductions$umap@cell.embeddings))
umap_df$combi <- factor(umap_df$combi, levels = c("protocol-driven - well-mixed",
                                                  "cell line-driven - well-mixed",
                                                  "cell line-driven - in vitro specific",
                                                  "protocol-driven - in vitro specific"))
set.seed(123)
p1 <- ggplot() + 
  geom_point(data = umap_df[sample(rownames(umap_df)),], aes(x = umap_1, y = umap_2, col = combi), alpha = 1, size = .3, stroke = 0) +
  scale_color_manual(values =  c("cell line-driven - in vitro specific" = "#1b9e76ff",
                                 "cell line-driven - well-mixed" = "#c18abbff",
                                 "protocol-driven - in vitro specific" = "#a2cd8dff",
                                 "protocol-driven - well-mixed" = "#6964aaff")) + 
  theme_void() +  
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=4), title="")) + NoLegend()

try(p1+ ggsave("figures/Fig3/Fig3E_UMAP_well_mixed_protocol_driven.png", 
             width = 3.5, 
             height = 3), silent = T)
                        
p <- ggplot() + 
  geom_point(data = umap_df, aes(x = umap_1, y = umap_2, col = combi), alpha = 1, size = .3, stroke = 0) +
  scale_color_manual(values = c("cell line-driven - in vivo specific" = "#1b9e76ff",
                                "cell line-driven - well-mixed" = "#c4afc8ff",
                                "protocol-driven - in vivo specific" = "#b9d5bbff",
                                "protocol-driven - well-mixed" = "#6964aaff")) + 
  theme_void() +  
  theme(text = element_text(size = 8, colour = "black"), 
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=4), title=""))  

pdf("figures/Fig3/Fig3E_UMAP_well_mixed_protocol_driven_legend.pdf", 
    width = 5, height = 3, pointsize = 8)

try(print(grid::grid.newpage() + 
            grid::grid.draw(cowplot::get_legend(p))), silent = T)
dev.off()

#### Fig 3F ####
df_new <- data.frame("division" = c(as.character(umap_invitro$annotation_coarse),
                                    as.character(umap_invitro$invivo_mixed)),
                     "protocol_driven" = c(umap_invitro$protocol_driven,
                                           umap_invitro$protocol_driven),
                     "in_vitro_specific" = c(umap_invitro$invivo_mixed,
                                             umap_invitro$invivo_mixed),
                     "rep" = c(rep(1,dim(umap_invitro)[1]),
                               rep(2,dim(umap_invitro)[1])),
                     "freq" = 1,
                     "ID" = c(rownames(umap_invitro),
                              rownames(umap_invitro)),
                     "annotation_coarse" = c(as.character(umap_invitro$annotation_coarse),
                                             as.character(umap_invitro$annotation_coarse)),
                     "annotation_fine" = c(umap_invitro$annotation_fine,
                                           umap_invitro$annotation_fine),
                     "protocol" = c(umap_invitro$protocol,
                                    umap_invitro$protocol)) 

df_new$annotation_coarse <- factor(df_new$annotation_coarse, 
                                   levels = names(annotation_coarse_colors))

df_new$division <- factor(df_new$division, 
                          levels = c(names(annotation_coarse_colors),
                                     c("well-mixed",
                                       "in vitro specific")))

df_new$protocol_driven <- as.character(df_new$protocol_driven)
df_new$protocol_driven[grepl("protocol_permissive",df_new$protocol_driven)] <- "protocol-driven"
df_new$protocol_driven[grepl("cell.line_biased",df_new$protocol_driven)] <- "cell line-driven"

df_new$protocol_driven <- factor(df_new$protocol_driven, 
                             levels = c("protocol-driven",
                                        "cell line-driven"))

p <- ggplot(df_new,
            aes(x = rep, 
                stratum = division, 
                alluvium = ID,
                y = freq,
                label = division, 
                fill = annotation_coarse))  + 
  scale_fill_manual(values = unlist(annotation_coarse_colors)) + 
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 4.1/6, alpha = 1) + 
  geom_stratum(width = 4.1/6, 
               col = "white", linewidth = 0) +
  geom_text(stat = "stratum", col = "white", size.unit = "pt", size = 6) + 
  scale_y_continuous(name = "number of cells (in thousands)",
                     labels = function(y) y / 1000) +
  xlab("") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  NoLegend() + 
  facet_wrap(~protocol_driven) 
try(p + ggsave("figures/Fig3/Fig3F_ALLUVIUM_well_mixed_protocol_driven.pdf", 
               width = 120, height = 90, units = "mm"), silent = T)
