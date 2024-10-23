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
dir.create("figures/SFig1", showWarnings = FALSE)

combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))
umap_df <- cbind(combined.obj@meta.data, combined.obj@reductions$umap@cell.embeddings)
expr_df <- cbind(t(combined.obj@assays$RNA$data), combined.obj@reductions$umap@cell.embeddings)

### Panel 1A ###
percellineandprotocol <- combined.obj@meta.data %>% group_by(protocol, cell.line.repl, replicate) %>%  tally()
percellineandprotocol$protocol <- factor(percellineandprotocol$protocol, rev(c("dorsal", "ventral", "striatum", "midbrain")))
percellineandprotocol$replicate <- factor(percellineandprotocol$replicate, levels = c("ORG1","ORG2"))
percellineandprotocol$cell.line_num <- as.numeric(percellineandprotocol$cell.line.repl) + 1/4*as.numeric(percellineandprotocol$replicate)

percellineandprotocol$cell.line_replicate <- paste0(percellineandprotocol$cell.line.repl,"_", percellineandprotocol$replicate)
percellineandprotocol$cell.line_replicate <- factor(percellineandprotocol$cell.line_replicate, levels = c("H1_ORG1","H1_ORG2",
                                                                                                          "H9_ORG1","H9_ORG2",
                                                                                                          "ROZH_5_ORG1","ROZH_5_ORG2",
                                                                                                          "XUJA_2_ORG1","XUJA_2_ORG2",
                                                                                                          "UOFV_1_ORG1","UOFV_1_ORG2",
                                                                                                          "176_E1_ORG1","176_E1_ORG2",
                                                                                                          "176_E2_ORG1","176_E2_ORG2",
                                                                                                          "177_ORG1","177_ORG2",
                                                                                                          "178_ORG1","178_ORG2"))
p <- ggplot(percellineandprotocol, aes(y = protocol, x = cell.line_num, label = round(n/1000,2), col = protocol)) + 
  geom_point(data = percellineandprotocol, aes(y = protocol, x = cell.line_num, size = round(n/1000,2), col = protocol)) + 
  scale_color_manual(values=combined.obj@misc$protocol_colors) + xlab("") + ylab("") + 
  theme_minimal() + 
  scale_x_continuous(limits = range(percellineandprotocol$cell.line_num), 
                     breaks = 1:length(levels(percellineandprotocol$cell.line.repl)) + .375, 
                     labels = levels(percellineandprotocol$cell.line.repl)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  geom_text(data = percellineandprotocol, vjust=-3, color = "black", cex = 1.5) + NoLegend()

p

try(p+ggsave("figures/SFig1/SuplFig1A_QC_cell_numbers.pdf", 
             width = 5, height = 2.5), silent = T)

### Panel 1B ###
marker_list <- 
  c("radial glia" = "PAX6",
    "radial glia" = "HOPX",
    "dorsal intermediate progenitors" = "EOMES", 
    
    "immature dorsal exc. neurons" = "SLA",
    "dorsal exc. neurons DL" = "TBR1",  
    "dorsal exc. neurons UL" = "SATB2",
    
    "ventral inh. neurons" = "DLX2", 
    "ventral inh. neurons" = "SST", 
    
    "floorplate prog." = "FOXA2", 
    "floorplate prog." = "EN1", 
    "medium spiny neurons" =  "SIX3", 
    "medium spiny neurons" =  "SP9", 
    
    "eye" = "RAX", 
    "fibroblasts"  = "COL3A1", 
    "muscle" = "DES" 
  )

featureplot_fun <- function(genename, gene_info) { 
  p <- ggplot(expr_df[order(expr_df[,genename]),
                      c("umap_1","umap_2",genename)], 
              aes_string(x = "umap_1", y = "umap_2", col = genename)) + 
    geom_point(size = 0.2, stroke = 0) + 
    theme_void() + 
    scale_colour_gradient(low = "snow2", high = "darkblue") +
    theme(text = element_text(size = 8, colour = "black"), 
          title = element_text(size = 8, colour = "black"),
          subtitle = element_text(size = 6, colour = "black"),
          legend.text = element_text(size = 8, colour = "black"),
          strip.text = element_text(size = 8, colour = "black"),
          legend.key.size = unit(.8, "lines")) + 
    labs(title=genename, subtitle=gene_info) +
    guides(col=guide_colorbar(title=""))  
  return(p)
}

p1 <- featureplot_fun(marker_list[1], names(marker_list)[1])
p2 <- featureplot_fun(marker_list[2], names(marker_list)[2])
p3 <- featureplot_fun(marker_list[3], names(marker_list)[3])
p4 <- featureplot_fun(marker_list[4], names(marker_list)[4])
p5 <- featureplot_fun(marker_list[5], names(marker_list)[5])
p6 <- featureplot_fun(marker_list[6], names(marker_list)[6])
p7 <- featureplot_fun(marker_list[7], names(marker_list)[7])
p8 <- featureplot_fun(marker_list[8], names(marker_list)[8])
p9 <- featureplot_fun(marker_list[9], names(marker_list)[9])
p10 <- featureplot_fun(marker_list[10], names(marker_list)[10])
p11 <- featureplot_fun(marker_list[11], names(marker_list)[11])
p12 <- featureplot_fun(marker_list[12], names(marker_list)[12])
p13 <- featureplot_fun(marker_list[13], names(marker_list)[13])
p14 <- featureplot_fun(marker_list[14], names(marker_list)[14])
p15 <- featureplot_fun(marker_list[15], names(marker_list)[15])

try(cowplot::plot_grid(p1,p2,p3,
                       p4,p5,p6,
                       p7,p8,p9,
                       p10,p11,p12,
                       p13,p14,p15,
                       ncol = 3) + 
      ggsave("figures/SFig1/SuplFig1B_UMAP_marker.png", width = 4.5, height = 7.5), silent = T)

### Panel 1C ###
p <- ggplot(umap_df, aes(fill = annotation_fine_num, x=annotation_coarse)) + 
  geom_bar() + ylab("number of cells") + 
  guides(fill=guide_legend(title="protocol")) + 
  scale_fill_manual(values = combined.obj@misc$annotation_fine_num_colors) +
  xlab(element_blank()) + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        text = element_text(size = 8, colour = "black"), 
        axis.text = element_text(size = 6, colour = "black"), 
        axis.title = element_text(size = 8, colour = "black"),
        title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8, colour = "black")) + 
  facet_wrap(~protocol) + NoLegend()
p

try(p + ggsave("figures/SFig1/SuplFig1C_barplot_celltypes_scaled.pdf", 
               width = 120, height = 72, units = "mm"), silent = T) 

### Panel 1D ### 
# voxhunt Analysis #

df <- structure_markers('E13', annotation_level = "custom_2")
regional_markers <- df%>%
  group_by(group) %>%
  top_n(15, auc) %>% 
  {unique(.$gene)}

vox_map <- voxel_map(combined.obj[,grepl("ventral|dorsal|midbrain", 
                                         combined.obj$annotation_coarse)], 
                     genes_use = regional_markers, 
                     group_name = "annotation_coarse", 
                     stage = "E13", 
                     method = "spearman")

pdf("figures/SFig1/SuplFig1D_voxhunt_by_celltype.pdf", width = 8, height = .9, pointsize = 8)
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
  (plot_map(vox_map, nrow = 1)[[5]] + theme(text = element_text(size = 8, colour = "black"), 
                                            title = element_text(size = 8, colour = "black"),
                                            legend.text = element_text(size = 8, colour = "black"),
                                            strip.text = element_text(size = 8, colour = "black"))) +
  (plot_map(vox_map, nrow = 1)[[6]] + theme(text = element_text(size = 8, colour = "black"), 
                                            title = element_text(size = 8, colour = "black"),
                                            legend.text = element_text(size = 8, colour = "black"),
                                            strip.text = element_text(size = 8, colour = "black"))) +
  patchwork::plot_layout(ncol = 6)
dev.off()

### Panel 1E ### 

p_H1_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                      aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + 
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "H1" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_H1_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                       aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() +# ggtitle("ventral") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "H1" & umap_df$protocol == "ventral",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_H1_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                        aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "H1" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_H1_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                        aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "H1"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

p_H9_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                      aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("H9") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "H9" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

# no cells:
p_H9_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                       aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void()

p_H9_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                        aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "H9" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_H9_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                        aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "H9"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

p_ROZH_5_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                          aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("ROZH_5") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "ROZH_5" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_ROZH_5_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                           aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("ventral") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "ROZH_5" & umap_df$protocol == "ventral",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_ROZH_5_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "ROZH_5" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_ROZH_5_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "ROZH_5"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

p_XUJA_2_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                          aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("XUJA_2") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "XUJA_2" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_XUJA_2_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                           aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("ventral") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "XUJA_2" & umap_df$protocol == "ventral",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_XUJA_2_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "XUJA_2" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_XUJA_2_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "XUJA_2"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

p_UOVF_1_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                          aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("UOFV_1") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "UOFV_1" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_UOVF_1_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                           aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("ventral") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "UOFV_1" & umap_df$protocol == "ventral",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_UOVF_1_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "UOFV_1" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_UOVF_1_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "UOFV_1"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

p_176_E1_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                          aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("176_E1") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "176_E1" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_176_E1_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                           aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("ventral") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "176_E1" & umap_df$protocol == "ventral",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_176_E1_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "176_E1" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_176_E1_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "176_E1"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

p_176_E2_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                          aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("176_E2") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "176_E2" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_176_E2_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                           aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("ventral") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "176_E2" & umap_df$protocol == "ventral",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_176_E2_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "176_E2" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_176_E2_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                            aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "176_E2"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

p_177_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                       aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("177") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "177" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_177_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                        aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("ventral") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "177" & umap_df$protocol == "ventral",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_177_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                         aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "177" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_177_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                         aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "177"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

p_178_dorsal <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                       aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("178") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "178" & umap_df$protocol == "dorsal",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_178_ventral <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                        aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("ventral") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "178" & umap_df$protocol == "ventral",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_178_midbrain <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                         aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("midbrain") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "178" & umap_df$protocol == "midbrain",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])
p_178_striatum <- ggplot(umap_df, vars = c("umap_1", "umap_2", "cell.line.repl"),
                         aes(x = umap_1, y = umap_2, col = protocol)) +
  theme_void() + #ggtitle("striatum") +
  geom_point(size = 0.1, alpha=1, col = "snow2") +
  geom_point(data = umap_df[umap_df$cell.line.repl == "178"  & umap_df$protocol == "striatum",], size = 0.1, alpha=1) + scale_color_manual(values = combined.obj@misc$protocol_colors)+ NoLegend()#, col = ggplotColours(9)[1])

try(cowplot::plot_grid(cowplot::plot_grid(p_UOVF_1_dorsal,
                                          p_UOVF_1_ventral,
                                          p_UOVF_1_midbrain,
                                          p_UOVF_1_striatum,nrow = 4),
                       cowplot::plot_grid(p_XUJA_2_dorsal,
                                          p_XUJA_2_ventral,
                                          p_XUJA_2_midbrain,
                                          p_XUJA_2_striatum,nrow = 4),
                       cowplot::plot_grid(p_H1_dorsal,
                                          p_H1_ventral,
                                          p_H1_midbrain,
                                          p_H1_striatum,nrow = 4),
                       cowplot::plot_grid(p_H9_dorsal,
                                          p_H9_ventral,
                                          p_H9_midbrain,
                                          p_H9_striatum,nrow = 4),
                       cowplot::plot_grid(p_ROZH_5_dorsal,
                                          p_ROZH_5_ventral,
                                          p_ROZH_5_midbrain,
                                          p_ROZH_5_striatum,nrow = 4),
                       cowplot::plot_grid(p_177_dorsal,
                                          p_177_ventral,
                                          p_177_midbrain,
                                          p_177_striatum,nrow = 4),
                       cowplot::plot_grid(p_178_dorsal,
                                          p_178_ventral,
                                          p_178_midbrain,
                                          p_178_striatum,nrow = 4),
                       cowplot::plot_grid(p_176_E1_dorsal,
                                          p_176_E1_ventral,
                                          p_176_E1_midbrain,
                                          p_176_E1_striatum,nrow = 4),
                       cowplot::plot_grid(p_176_E2_dorsal,
                                          p_176_E2_ventral,
                                          p_176_E2_midbrain,
                                          p_176_E2_striatum,nrow = 4),
                       nrow = 1) + ggsave("figures/SFig1/SuplFig1E_UMAP_by_protocol_cell.line.png", width = 12, height = 4), silent = T)
