rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(limma)
library(edgeR)
library(xtable)
library(ggplot2)
library(topGO)
library(ggrepel)
library(dplyr)
library(Seurat)

source("./helper_functions/bulk_functions.R")
source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./bulkRNAseq_analysis/data/"
dir.create("figures/Fig5", showWarnings = FALSE)
dir.create("figures/SFig5", showWarnings = FALSE)

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

# get transcription factors
library(biomaRt)
library(org.Hs.eg.db)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                   host="asia.ensembl.org")

tf_list <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                 filters = 'go',
                 values = 'GO:0003700', # GO term for transcription factor activity
                 mart = ensembl)

#### DGE ANALYIS - all time points ####
sample_sele = metadata$protocol != "pluripotent"

d0 = DGEList(counts[,sample_sele])
cutoff <- 20
drop <- which(apply(cpm(d0), 1, function(x) {sum(x >= 5)}) < cutoff)
d <- d0[-drop,] 

protocol = metadata$protocol[sample_sele]
protocol = factor(protocol, levels = levels(protocol)[which(levels(protocol) %in% protocol)])

time_num = metadata$time_num[sample_sele]
time_num = factor(time_num, levels = sort(unique(time_num)))

cell.line = metadata$cell.line[sample_sele]
cell.line = factor(cell.line, levels = levels(cell.line)[which(levels(cell.line) %in% cell.line)])

mm = model.matrix(~ 0 + protocol + time_num + cell.line)

colnames(mm) = gsub(':', '.', colnames(mm))
y = voom(d, mm, plot=TRUE)
fit = lmFit(y, mm)

genes_list = vector()
protocol_list = vector()
adj_p_value_list = vector()
log_fc_list = vector()
number_of_genes_plot = 1

for(protocol_i in names(protocol_colors)[-5]){
  
  # test each protocol against the average of the other 3
  if(protocol_i == 'dorsal') {
    contr = makeContrasts(
      MvsAll = protocoldorsal - (protocolventral + protocolmidbrain + protocolstriatum)/3,
      levels=colnames(mm)
    )
  } else if(protocol_i == 'ventral') {
    contr = makeContrasts(
      MvsAll = protocolventral - (protocoldorsal + protocolmidbrain + protocolstriatum)/3,
      levels=colnames(mm)
    )
  } else if(protocol_i == 'midbrain') {
    contr = makeContrasts(
      MvsAll = protocolmidbrain - (protocolventral + protocoldorsal + protocolstriatum)/3,
      levels=colnames(mm)
    )
  } else if(protocol_i == 'striatum') {
    contr = makeContrasts(
      MvsAll = protocolstriatum - (protocolventral + protocolmidbrain + protocoldorsal)/3,
      levels=colnames(mm)
    )
  } else next
  
  tmp = contrasts.fit(fit, contr)
  tmp = eBayes(tmp)
  
  tt = topTable(tmp, sort.by='p', p.value=0.01, adjust.method = "BH", number = 2000)
  tt = tt[tt$logFC>0.2,]
  
  if(1) {
    df_volcano = data.frame(p_value=-log10(tmp$p.value[,1]), logfc=tmp$coefficients[,1],
                            names=map_symbols(rownames(tmp), gene_symbols))
    
    sorted_df = df_volcano[order(df_volcano$p_value, decreasing=TRUE),]
    top_genes = sorted_df$names[sorted_df$logfc>0.2][1:10]
    
    p1 <- ggplot(df_volcano, aes(x=logfc, y=p_value, label=names)) +
      geom_point(color=ifelse(df_volcano$names %in% top_genes,'red','black')) +
      ggrepel::geom_label_repel(aes(label=ifelse(names %in% top_genes,as.character(names),'')),
                                color='red', 
                                min.segment.length = .3, 
                                label.size = NA, 
                                fill = "transparent", 
                                box.padding = .3, 
                                label.padding = 0, 
                                max.overlaps = Inf) +
      xlim(min(df_volcano$logfc)-.3, max(df_volcano$logfc)+0.5) +
      xlab('log2FoldChange') + ylab('-log10(p-value)') +
      theme_minimal() +
      theme(text=element_text(size=20))
    print(p1)
    
    try(p1 + ggsave(paste0('figures/Fig5/Fig5A_', protocol_i, '_volcano.pdf'), 
                    width = 4.5, height = 4, dpi = 600), silent = T)
  }
  
  # create plot of expression over time for a gene in each protocol and cell line
  for(gene in rownames(tt)[1:number_of_genes_plot]) {
    
    if(protocol_i == "striatum") { 
      gene = rownames(tt)[2]
    }
    expression = expr[gene,]
    plot_df = data.frame(expression = as.numeric(expression), 
                         metadata)
    
    # average over replicates
    ag = aggregate(plot_df$expression, 
                   by=list(plot_df$protocol, 
                           plot_df$time_num, 
                           plot_df$cell.line), mean)
    
    colnames(ag) <- c("protocol", "time_num","cell.line","expression")
    
    # use plurip data as timepoints 0 in each protocol
    plurip = ag[ag$protocol=='pluripotent',]
    plurip$protocol = rep('dorsal', length(plurip$protocol))
    ag = rbind(ag, plurip)
    plurip$protocol = rep('ventral', length(plurip$protocol))
    ag = rbind(ag, plurip)
    plurip$protocol = rep('midbrain', length(plurip$protocol))
    ag = rbind(ag, plurip)
    plurip$protocol = rep('striatum', length(plurip$protocol))
    ag = rbind(ag, plurip)
    ag = ag[ag$protocol!='pluripotent',]
    
    symbol = map_symbols(ensembl=gene,map=gene_symbols)
    
    ag$protocol = factor(ag$protocol, levels=c('dorsal', 'ventral', 'midbrain', 'striatum'))
    
    ag$time_str <- as.character(ag$time_num)
    ag$time_str <- sub("13","13/16", ag$time_str)

    p <- ggplot(ag, aes(x=time_num, y=expression, color=cell.line)) +
      geom_line() + scale_color_manual(values = cell.line_colors) +
      facet_wrap(~protocol) +
      ggtitle(symbol) +
      geom_point(size=4) +
      theme_minimal() +
      scale_x_continuous(limits = range(ag$time_num),
                         breaks = unique(sort(ag$time_num)),
                         labels = unique(ag$time_str[order(ag$time_num)])) +
      theme(text=element_text(size=20),
            axis.text.x = element_text(angle = 45,hjust=1)) + xlab("time (days)")
    print(p)
      
    filename = paste0('figures/Fig5/Fig5B_', protocol_i, '_', symbol, '_', gene, '.pdf')
    try(p+ggsave(filename, width = 8, height = 5), silent = T)
  }
  
  if(1) {
    df_volcano = data.frame(p_value=-log10(tmp$p.value[,1]), logfc=tmp$coefficients[,1],
                            names=map_symbols(rownames(tmp), gene_symbols))
    
    tf_num = ifelse(df_volcano$names %in% tf_list$hgnc_symbol, 1, 0)
    df_volcano = df_volcano[order(tf_num),]
    
    sorted_df = df_volcano[order(df_volcano$p_value, decreasing=TRUE),]
    top_genes = sorted_df$names[(sorted_df$logfc>0.2) & (sorted_df$names %in% tf_list$hgnc_symbol)][1:10]
    
    is_tf = df_volcano$names %in% tf_list$hgnc_symbol
    is_top = df_volcano$names %in% top_genes
    color_vector = rep('grey', nrow(df_volcano))
    color_vector[is_tf] = 'black'
    color_vector[is_top] = 'red'
    p1 <- ggplot(df_volcano, aes(x=logfc, y=p_value, label=names)) +
      geom_point(color=color_vector, alpha=ifelse(is_tf, 1, .5)) +
      ggrepel::geom_label_repel(aes(label=ifelse(is_top, as.character(names),'')),
                                color='red', 
                                min.segment.length = .3, 
                                label.size = NA, 
                                fill = "transparent", 
                                box.padding = .3, 
                                label.padding = 0, 
                                max.overlaps = Inf) +
      xlim(min(df_volcano$logfc)-.3, max(df_volcano$logfc)+0.5) +
      xlab('log2FoldChange') + ylab('-log10(p-value)') +
      theme_minimal() 
    
    print(p1)
    try(p1 + ggsave(paste0('figures/SFig5/SuplFig5B_', protocol_i, '_tf_volcano.pdf'), 
                    width = 4.5, height = 4, dpi = 600), silent = T)
  }
  
  adj_p_value_list = append(adj_p_value_list, tt$adj.P.Val)
  genes_list = append(genes_list, map_symbols(rownames(tt), gene_symbols))
  protocol_list = append(protocol_list, rep(protocol_i, length(rownames(tt))))
  log_fc_list = append(log_fc_list, tt$logFC)
}

df_alltimes = data.frame(Gene = genes_list, 
                adj_p_value = adj_p_value_list, 
                protocol = protocol_list,
                logFC = log_fc_list)

write.csv(df_alltimes, 
          paste0(DataDir,'bulk_RNAseq_marker.csv'))

### DGE ANALYIS - day 40 - protocol permissive cell lines ####

genes_list = vector()
protocol_list = vector()
adj_p_value_list = vector()
log_fc_list = vector()
number_of_genes_plot = 1

for(protocol_i in names(protocol_colors)[-5]){
  
  sample_sele = (metadata$protocol != "pluripotent" & 
                   metadata$time_num == 40)
  
  d0 = DGEList(counts[,sample_sele])
  cutoff <- 10 # less samples for single time points
  drop <- which(apply(cpm(d0), 1, function(x) {sum(x >= 5)}) < cutoff)
  d <- d0[-drop,] 
  
  protocol = metadata$protocol[sample_sele]
  protocol = factor(protocol, levels = levels(protocol)[which(levels(protocol) %in% protocol)])
  
  cell.line = metadata$cell.line[sample_sele]
  cell.line = factor(cell.line, levels = levels(cell.line)[which(levels(cell.line) %in% cell.line)])
  
  growth = metadata$growth[sample_sele] 
  growth[growth == "protocol-driven"] = "unbiased"
  growth[growth == "cell line-driven"] = "biased"
  growth = as.factor(growth)
  
  protocol_growth = paste0(protocol,"_",growth)
  protocol_growth = as.factor(protocol_growth)  
  
  mm = model.matrix(~ 0 + protocol_growth + cell.line)
  
  colnames(mm) = gsub(':', '.', colnames(mm))
  y = voom(d, mm, plot=TRUE)
  fit = lmFit(y, mm)
  
  # test each protocol against the average of the other 3
  if(protocol_i == 'dorsal') {
    contr = makeContrasts(
      MvsAll = protocol_growthdorsal_unbiased - 
        (protocol_growthdorsal_biased +  
           protocol_growthventral_unbiased +
           protocol_growthventral_biased +
           protocol_growthmidbrain_unbiased +
           protocol_growthstriatum_unbiased +
           protocol_growthstriatum_biased)/6,
      levels=colnames(mm)
    )
  } else if(protocol_i == 'ventral') {
    contr = makeContrasts(
      MvsAll = protocol_growthventral_unbiased - 
        (protocol_growthventral_biased +  
           protocol_growthdorsal_unbiased +
           protocol_growthdorsal_biased +
           protocol_growthmidbrain_unbiased +
           protocol_growthstriatum_unbiased +
           protocol_growthstriatum_biased)/6,
      levels=colnames(mm)
    )
  } else if(protocol_i == 'midbrain') {
    contr = makeContrasts(
      MvsAll = protocol_growthmidbrain_unbiased - 
        (protocol_growthventral_unbiased +  
           protocol_growthventral_biased +  
           protocol_growthdorsal_unbiased +
           protocol_growthdorsal_biased +
           protocol_growthstriatum_unbiased +
           protocol_growthstriatum_biased)/6,
      levels=colnames(mm)
    )
  } else if(protocol_i == 'striatum') {
    contr = makeContrasts(
      MvsAll = protocol_growthstriatum_unbiased - 
        (protocol_growthventral_unbiased +  
           protocol_growthventral_biased +  
           protocol_growthdorsal_unbiased +
           protocol_growthdorsal_biased +
           protocol_growthmidbrain_unbiased +
           protocol_growthstriatum_biased)/6,
      levels=colnames(mm)
    )
  } else next
  
  tmp = contrasts.fit(fit, contr)
  tmp = eBayes(tmp)
  
  tt = topTable(tmp, sort.by='p', p.value=0.01, adjust.method = "BH", number = 2000)
  tt = tt[tt$logFC>0.2,]
  
  if(1) {
    df_volcano = data.frame(p_value=-log10(tmp$p.value[,1]), logfc=tmp$coefficients[,1],
                            names=map_symbols(rownames(tmp), gene_symbols))
    
    sorted_df = df_volcano[order(df_volcano$p_value, decreasing=TRUE),]
    top_genes = sorted_df$names[sorted_df$logfc>0.2][1:10]
    
    p1 <- ggplot(df_volcano, aes(x=logfc, y=p_value, label=names)) +
      geom_point(color=ifelse(df_volcano$names %in% top_genes,'red','black')) +
      ggrepel::geom_label_repel(aes(label=ifelse(names %in% top_genes,as.character(names),'')),
                                color='red', 
                                min.segment.length = .3, 
                                label.size = NA, 
                                fill = "transparent", 
                                box.padding = .3, 
                                label.padding = 0, 
                                max.overlaps = Inf) +
      xlim(min(df_volcano$logfc)-.3, max(df_volcano$logfc)+0.5) +
      xlab('log2FoldChange') + ylab('-log10(p-value)') +
      theme_minimal() +
      theme(text=element_text(size=20))
    print(p1)
    
    try(p1 + ggsave(paste0('figures/SFig5/SuplFig5B_', protocol_i, '_volcano_day40_permissive.pdf'), 
                    width = 4.5, height = 4, dpi = 600), silent = T)
  }
  
  adj_p_value_list = append(adj_p_value_list, tt$adj.P.Val)
  genes_list = append(genes_list, map_symbols(rownames(tt), gene_symbols))
  protocol_list = append(protocol_list, rep(protocol_i, length(rownames(tt))))
  log_fc_list = append(log_fc_list, tt$logFC)
}

df_day40_protocoldriven = data.frame(Gene = genes_list, 
                adj_p_value = adj_p_value_list, 
                protocol = protocol_list,
                logFC = log_fc_list)

df_list <- list("marker_alltimes" = df_alltimes, 
                "marker_day40_protocoldriven" = df_day40_protocoldriven)

writexl::write_xlsx(df_list, 
                    "tables/SuplTable3_DEG_bulkRNA_invitro.xlsx")
