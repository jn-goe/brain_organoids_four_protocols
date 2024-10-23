rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(limma)
library(edgeR)
library(xtable)
library(ggplot2)
library(topGO)
library(ggrepel)
library(dplyr)
library(Seurat)
library(DESeq2)
library(pheatmap)

source("./helper_functions/colors_metadata.R")
source("./helper_functions/bulk_functions.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./bulkRNAseq_analysis/data/"
dir.create("figures/Fig4", showWarnings = FALSE)
dir.create("figures/SFig4", showWarnings = FALSE)

cell.line_colors <- cell.line_colors[-1]

#### DATA LOADING AND PREPROCESSING
gene_symbols = read.table(paste0(DataDir,"ii.GRCh38_202211_mapTable.tsv"), header=TRUE)
colnames(gene_symbols) = c('ensembl', 'symbol')

metadata <- read.table(paste0(DataDir,"bulkRNAseq_metadata.csv"), sep=',', header=TRUE)
expr <- readRDS(paste0(DataDir,"allsamples_counts_TPM.RDS"))
counts <- readRDS(paste0(DataDir,"allsamples_counts.RDS"))

# EXCLUDE SAMPLES WITH LESS THAN 100,000 READS
metadata$depth = apply(counts, 2, sum)
expr = expr[,metadata$depth > 100000]
counts = counts[,metadata$depth > 100000]
metadata = metadata[metadata$depth > 100000,]

metadata$cell.line <- factor(metadata$cell.line, levels = names(cell.line_colors))
metadata$protocol <- factor(metadata$protocol, levels = names(protocol_colors))

#### COMPUTE VARIABLE GENES ####
genes_all = map_symbols(rownames(expr), map=gene_symbols)

vsd = vst(as.matrix(counts))
outlier_list = c(genes_all[grepl("^RPL|^RS|^RPS|^MRP", genes_all)], 
                 genes_all[grepl("^MT\\.|^MT-", genes_all)])
vsd = vsd[!genes_all %in% outlier_list,]

n_genes = 2000

for(timepoint in unique(metadata$time_str)) { 
  
  sele = metadata$time_str == timepoint
  expr_sele = vsd[,sele]
  meta_sele = metadata[sele,]
  
  var = apply(expr_sele, 1, var)
  
  var_genes = expr_sele[order(var, decreasing=TRUE)[1:n_genes],]
  rownames(var_genes) = map_symbols(rownames(var_genes), map=gene_symbols)
  rownames(meta_sele) = colnames(var_genes)
  
  anno_df = data.frame(growth=meta_sele$growth,
                       cell.line = meta_sele$cell.line,
                       protocol = meta_sele$protocol)
  rownames(anno_df) = colnames(var_genes)
  
  ann_colors = list(protocol=protocol_colors, 
                    cell.line=cell.line_colors, 
                    growth=growth_colors)
  
  mat <- as.matrix(var_genes)
  rownames(mat) <- make.unique(rownames(mat))
  
  res <- pheatmap(mat,
                  scale = 'row',
                  show_colnames = F,
                  show_rownames = F,
                  treeheight_row = 0)
  dev.off()
  
  mat_new <- mat[res$tree_row$order, ]
  
  if(timepoint == "day000") { 
    
    p <- "pluripotent"
    anno_df_in <- anno_df[anno_df$protocol == p,]

    mat_sub <- mat_new
    mat_sub <- mat_sub[,order(anno_df_in$cell.line)]
    anno_df_in <- anno_df_in[order(anno_df_in$cell.line),]
    
    mat_sub_new <- matrix(NA, nrow = dim(mat_sub)[1], ncol = length(levels(metadata$cell.line)))
    rownames(mat_sub_new) <- rownames(mat_sub)
    colnames(mat_sub_new) <- levels(metadata$cell.line)
    
    for(c_in in unique(anno_df_in$cell.line)) { 
      mat_sub_new[,c_in] <- rowMeans(mat_sub[,anno_df_in$cell.line == c_in, drop = F])
    }
    
    in_mat <- cor(mat_sub_new)
    
    if(any(is.na(in_mat))) { 
      ind <- which(table(which(is.na(in_mat), arr.ind = T)[,1])>2)
      for(i in ind) { 
        in_mat[ind,ind] <- NA
      }
    }
    colnames(in_mat) <- rownames(in_mat) <- c(colnames(mat_sub_new))#,"all")
    
    ann_new <- data.frame("cell.line" = levels(metadata$cell.line))
    rownames(ann_new) <- levels(metadata$cell.line)
    
    ann_colors <- list("cell.line" = cell.line_colors)
    
    pdf(paste0("figures/Fig4/Fig4H_time_",timepoint,"_",p,"_cor_heatmap.pdf"), 
        width = 3.6, height = 3)
    pheatmap(in_mat, 
             colorRampPalette(c("grey80", "grey20"))(100), na_col = "transparent",
             cluster_cols = F, 
             cluster_rows = F,
             main = paste0(timepoint," - ",p, " protocol")
    )
    dev.off()    
    
    anno_df <- anno_df_in[,c("cell.line"), drop = F]

    pdf(paste0("figures/Fig4/Fig4F_time_",timepoint,"_expression_heatmap.pdf"), width = 10, height = 6)
    pheatmap(mat_new,
             annotation_col = anno_df,
             scale = 'row',
             color = pals::coolwarm(50),
             annotation_colors = ann_colors,
             treeheight_row = 0,
             show_colnames = F,
             show_rownames = F,
             cluster_rows = F,
             fontsize=18, 
             main=timepoint)
    dev.off()
    
    # saveRDS(list(mat_new, anno_df, ann_colors), 
    #         file = paste0(DataDir,timepoint,"_expressionmap.RDS"))
  } else { 
    for(p in unique(anno_df$protocol)) { 
      
      anno_df_in <- anno_df[anno_df$protocol == p,]

      mat_sub <- mat_new[,anno_df$protocol == p]
      mat_sub <- mat_sub[,order(anno_df_in$cell.line)]
      anno_df_in <- anno_df_in[order(anno_df_in$cell.line),]
      
      mat_sub_new <- matrix(NA, nrow = dim(mat_sub)[1], ncol = length(levels(metadata$cell.line)))
      rownames(mat_sub_new) <- rownames(mat_sub)
      colnames(mat_sub_new) <- levels(metadata$cell.line)
      
      for(c_in in unique(anno_df_in$cell.line)) { 
        mat_sub_new[,c_in] <- rowMeans(mat_sub[,anno_df_in$cell.line == c_in, drop = F])
      }
      
      in_mat <- cor(mat_sub_new)
      
      if(any(is.na(in_mat))) { 
        ind <- which(table(which(is.na(in_mat), arr.ind = T)[,1])>2)
        for(i in ind) { 
          in_mat[ind,ind] <- NA
        }
      }
      colnames(in_mat) <- rownames(in_mat) <- c(colnames(mat_sub_new))#,"all")
      
      ann_new <- data.frame("cell.line" = levels(metadata$cell.line))
      rownames(ann_new) <- levels(metadata$cell.line)

      if(timepoint %in% c("day040","day120") & p == "dorsal") { 
        pdf(paste0("figures/Fig4/Fig4I_time_",timepoint,"_protocol_",p,"_cor_heatmap.pdf"), 
            width = 3.6, height = 3)
        pheatmap(in_mat,
                 colorRampPalette(c("grey80", protocol_colors[[p]]))(100), na_col = "transparent",
                 cluster_cols = F, 
                 cluster_rows = F,
                 main = paste0(timepoint," - ",p, " protocol")
        )
        dev.off()             
      } else {
        if(p == "dorsal") {
          pdf(paste0("figures/SFig4/SuplFig4D_time_",timepoint,"_protocol_",p,"_cor_heatmap.pdf"), 
              width = 3.6, height = 3)
        } else if(p == "ventral") {
          pdf(paste0("figures/SFig4/SuplFig4E_time_",timepoint,"_protocol_",p,"_cor_heatmap.pdf"), 
              width = 3.6, height = 3)
        } else if(p == "midbrain") {
          pdf(paste0("figures/SFig4/SuplFig4F_time_",timepoint,"_protocol_",p,"_cor_heatmap.pdf"), 
              width = 3.6, height = 3)
        } else if(p == "striatum") {
          pdf(paste0("figures/SFig4/SuplFig4G_time_",timepoint,"_protocol_",p,"_cor_heatmap.pdf"), 
              width = 3.6, height = 3)
        }
        pheatmap(in_mat,
                 colorRampPalette(c("grey80", protocol_colors[[p]]))(100), na_col = "transparent",
                 cluster_cols = F, 
                 cluster_rows = F,
                 main = paste0(timepoint," - ",p, " protocol")
        )
        dev.off()             
      }  
    }
    
    pdf(paste0("figures/Fig4/Fig4G_time_",timepoint,"_heatmap.pdf"), width = 10, height = 6)
    pheatmap(mat_new,
             annotation_col = anno_df,
             scale = 'row',
             color = pals::coolwarm(50),
             annotation_colors = ann_colors,
             treeheight_row = 0,
             show_colnames = F,
             show_rownames = F,
             cluster_rows = F,
             fontsize=15,
             main=timepoint)
    dev.off()

  }
}
