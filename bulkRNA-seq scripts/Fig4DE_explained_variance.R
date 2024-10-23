rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(dplyr)
library(ggplot2)
library(Seurat)
library(variancePartition)
library(edgeR)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
setwd("./")
DataDir <- "./bulkRNAseq_analysis/data/"
dir.create("figures/Fig4", showWarnings = FALSE)

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

### Fig4D ####
#---EXPL VARIANCE PIEPLOT ---#
sample_sele = metadata$protocol != "pluripotent"

d0 = DGEList(counts[,sample_sele])
cutoff <- 20
drop <- which(apply(cpm(d0), 1, function(x) {sum(x >= 5)}) < cutoff)
d <- d0[-drop,] 

protocol = metadata$protocol[sample_sele]
time = factor(metadata$time_num[sample_sele], 
                  levels = sort(unique(metadata$time_num[sample_sele])))
cell.line = metadata$cell.line[sample_sele]

mm = model.matrix(~ 0 + protocol + time + cell.line)

info = data.frame(protocol, cell.line, time)
form = ~ (1|protocol) + (1|cell.line) + (1|time)
mm = model.matrix(~ cell.line + protocol + time)
colnames(mm) = gsub(':', '.', colnames(mm))

y = voom(d, mm, plot=TRUE)

expr_sele = expr[-drop, sample_sele]
abs_var = apply(expr_sele, 1, var)

varPart = fitExtractVarPartModel(y, form, info)
apply(varPart, 2, mean)

abs_var_part = varPart * abs_var
weighted_var_part = apply(abs_var_part, 2, sum) / sum(abs_var_part)

df = data.frame(factor=names(weighted_var_part),
                var=weighted_var_part,
                dummy=rep('1', length(weighted_var_part)),
                label=paste0(sprintf("%.0f", weighted_var_part*100),"%"))

df=df[order(df$var, decreasing = TRUE),]
prop = df$var/sum(df$var)

round(df$var/sum(df$var)*100,2)
#[1] 39.62 39.04 13.29  8.05

ypos = cumsum(prop) - 0.5*prop

p_pie <- ggplot(data=df, aes(fill=factor, y=var, x="")) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(y=ypos, label=label), size=5) +
  scale_fill_manual(values=condition_colors)
p_pie
try(p_pie + ggsave("figures/Fig4/Fig4D_expl_var_pie.pdf", width = 4, height = 4), silent = T)

### Fig4E ####
#---EXPL VARIANCE OVER TIME ---#

expl_var_time = vector()
var_vector = vector()
time_vector = vector()

times = metadata$time_num
protocol = metadata$protocol
cell.line = metadata$cell.line

d0 = DGEList(counts)
cutoff <- 20
drop <- which(apply(cpm(d0), 1, function(x) {sum(x >= 5)}) < cutoff)
d <- d0[-drop,] 
dim(d)

for(timepoint in unique(times)) {
  
  sample_time = times == timepoint
  
  if(timepoint != 0) { 
    protocol_time = protocol[sample_time]
    cell.line_time = cell.line[sample_time]
    
    info = data.frame(protocol_time, cell.line_time)
    form = ~ (1|protocol_time) + (1|cell.line_time)
    mm = model.matrix(~ protocol_time + cell.line_time)
  } else {
    cell.line_time = cell.line[sample_time]
    
    info = data.frame(cell.line_time)
    form = ~ (1|cell.line_time)
    mm = model.matrix(~ cell.line_time)
  }
  
  colnames(mm) = gsub(':', '.', colnames(mm))
  y = voom(d[,sample_time], mm, plot=TRUE)
  abs_var = apply(expr[-drop,sample_time], 1, var)

  varPart = fitExtractVarPartModel(y, form, info)
  abs_var_part = varPart * abs_var
  weighted_var_part = apply(abs_var_part, 2, sum) / sum(abs_var_part)
  
  if(timepoint == 0) {
    weighted_var_part['protocol_time'] <- 0
  }
  
  expl_var_time = append(expl_var_time, weighted_var_part['protocol_time'])
  expl_var_time = append(expl_var_time, weighted_var_part['cell.line_time'])
  expl_var_time = append(expl_var_time, weighted_var_part['Residuals'])
  
  var_vector = append(var_vector, 'protocol')
  var_vector = append(var_vector, 'cell.line')
  var_vector = append(var_vector, 'Residuals')
  
  time_vector = append(time_vector, rep(timepoint, 3))
}

df = data.frame(variable=var_vector,
                time=time_vector,
                expl_variance=expl_var_time)

time_25 = df[df$time==25,]
df = df[df$time!=25,]

df$variable <- factor(df$variable, levels = names(condition_colors))

df$time_str <- as.character(df$time)
df$time_str[grepl("13", df$time_str)] <- "13/16"

p_time <- ggplot(df, aes(x=time, y=expl_variance*100)) +
  geom_area(aes(fill=variable)) + 
  geom_point(data=time_25, aes(x=time, y=expl_variance*100,
                               fill=variable), show.legend=FALSE,
             shape=23, color="black", size=3) +
  theme_minimal() +
  scale_x_continuous(limits = range(df$time),
                     breaks = c(df$time,25),
                     labels = c(df$time_str,25)) +
  scale_fill_manual(values=condition_colors) +
  theme(text = element_text(size=20)) +
  ylab("explained variance [%]") + xlab("time (days)")
p_time
try(p_time + ggsave("figures/Fig4/Fig4E_expl_var_over_time.pdf", width = 7, height = 4), silent = T)
