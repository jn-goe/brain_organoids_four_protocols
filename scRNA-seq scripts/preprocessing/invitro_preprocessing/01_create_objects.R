rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T) 
library(Seurat)

# Setup ----------------------------
OutDirOrig <- "./scRNAseq_analysis/separate_organoids/cellranger/"
setwd(OutDirOrig)

allnames  <- list("31" = c("111-H1-D-ORG1","211-H9-D-ORG1","311-ROZH_5-D-ORG1","411-XUJA_2-D-ORG1"),
                  "32" = c("121-H1-V-ORG1","221-H9-V-ORG1","321-ROZH_5-V-ORG1","421-XUJA_2-V-ORG1"),
                  "33" = c("131-H1-M-ORG1","231-H9-M-ORG1","331-ROZH_5-M-ORG1","431-XUJA_2-M-ORG1"),
                  "34" = c("141-H1-S-ORG1","241-H9-S-ORG1","341-ROZH_5-S-ORG1","441-XUJA_2-S-ORG1"),
                  "35" = c("112-H1-D-ORG2","212-H9-D-ORG2","312-ROZH_5-D-ORG2","412-XUJA_2-D-ORG2"),
                  "36" = c("122-H1-V-ORG2","222-H9-V-ORG2","322-ROZH_5-V-ORG2","422-XUJA_2-V-ORG2"),
                  "37" = c("132-H1-M-ORG2","232-H9-M-ORG2","332-ROZH_5-M-ORG2","432-XUJA_2-M-ORG2"),
                  "38" = c("142-H1-S-ORG2","242-H9-S-ORG2","342-ROZH_5-S-ORG2","442-XUJA_2-S-ORG2"),
                  "41" = c("511-UOFV_1-D-ORG1","512-UOFV_1-D-ORG2","611-176_E1-D-ORG1","612-176_E1-D-ORG2"),
                  "42" = c("521-UOFV_1-V-ORG1","522-UOFV_1-V-ORG2","621-176_E1-V-ORG1","622-176_E1-V-ORG2"),
                  "43" = c("531-UOFV_1-M-ORG1","532-UOFV_1-M-ORG2","631-176_E1-M-ORG1","632-176_E1-M-ORG2"),
                  "44" = c("541-UOFV_1-S-ORG1","542-UOFV_1-S-ORG2","641-176_E1-S-ORG1","642-176_E1-S-ORG2"),
                  "69" = c("721-176_E2-V-ORG1","722-176_E2-V-ORG2","821-177-V-ORG1","822-177-V-ORG2","921-178-V-ORG1","922-178-V-ORG2"),
                  "70" = c("731-176_E2-M-ORG1","732-176_E2-M-ORG2","831-177-M-ORG1","832-177-M-ORG2","931-178-M-ORG1","932-178-M-ORG2"),
                  "71" = c("741-176_E2-S-ORG1","742-176_E2-S-ORG2","841-177-S-ORG1","842-177-S-ORG2","941-178-S-ORG1","942-178-S-ORG2"),
                  "76" = c("711-176_E2-D-ORG1","712-176_E2-D-ORG2","811-177-D-ORG1","911-178-D-ORG1","912-178-D-ORG2"))
cell.line <- c("H1","H9","ROZH_5","XUJA_2","UOFV","176","177","178") 

folderlist <- list.files(OutDirOrig)
folderlist <- folderlist[!grepl(".Rds|.RDS", folderlist)]
folderlist <- folderlist[!grepl(".csv", folderlist)]

for(lane in folderlist) {
  
  obj <- Read10X(data.dir = paste0(OutDirOrig,lane,"/filtered_feature_bc_matrix"))
  orig.idents <- allnames[[lane]]
  folderlist.sub <- list.files(paste0(OutDirOrig,lane))
  
  if("HTO_demux.csv" %in% folderlist.sub) {
    HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(obj$`Gene Expression`), 
                                                 unmapped.as.na = FALSE, 
                                                 map = NULL, 
                                                 species = "human")
    rownames(obj$`Gene Expression`) <- make.unique(HGNC.updated$Suggested.Symbol)
    
    obj <- CreateSeuratObject(obj$`Gene Expression`, 
                              min.cells = 0, 
                              min.features = 0)
    obj$batch <- lane
    
    demultiplex <- read.csv(paste0(OutDirOrig,lane,"/HTO_demux.csv"), sep = ",")
    demultiplex <- subset(demultiplex, hto_demux != "Doublet" & hto_demux != "Negative")
    demultiplex$hto_demux <- sub("HTO-", "", demultiplex$hto_demux)
    demultiplex$hto_demux <- sub(paste0("sSL00",lane,"_"), "", demultiplex$hto_demux)
    
    demultiplex$index <- paste0(demultiplex$index, "-1")
    demultiplex <- demultiplex[toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"-Si")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"-Di")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"-Vi")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"-Mi")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"i")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"-Sii")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"-Dii")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"-Vii")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"-Mii")) |
                                 toupper(demultiplex$hto_demux) %in% toupper(paste0(cell.line,"ii")),]
    obj <- subset(obj, cells = demultiplex$index)
    cellnames <- sub("-Sii", "", demultiplex$hto_demux)
    cellnames <- sub("-Dii", "", cellnames)
    cellnames <- sub("-Vii", "", cellnames)
    cellnames <- sub("-Mii", "", cellnames)
    cellnames <- sub("ii", "", cellnames)
    cellnames <- sub("-Si", "", cellnames)
    cellnames <- sub("-Di", "", cellnames)
    cellnames <- sub("-Vi", "", cellnames)
    cellnames <- sub("-Mi", "", cellnames)
    cellnames <- sub("i", "", cellnames)
    
    obj$cellline <- toupper(cellnames)
    obj$replicate <- "ORG1"
    obj$replicate[demultiplex$index[grepl("ii", demultiplex$hto_demux)]] <- "ORG2"
    
    Idents(obj) <- obj$cellline
    ls.Seurat <- SplitObject(obj)
    
    for(s in ls.Seurat) {
      cellline.ident <- orig.idents[grepl(unique(s$cellline),orig.idents)]
      s$orig.ident <- cellline.ident[grepl("ORG1",cellline.ident)]
      s$orig.ident[s$replicate == "ORG2"] <- cellline.ident[grepl("ORG2",cellline.ident)]
      
      if(unique(substr(s$orig.ident, 1,1)) == 5) {
        s$cellline <- "UOFV_1"
      }
      if(unique(substr(s$orig.ident, 1,1)) == 6) {
        s$cellline <- "176_E1"
      }
      if(unique(substr(s$orig.ident, 1,1)) == 7) {
        s$cellline <- "176_E2"
      }
      
      substring <- substr(s$orig.ident,5,length(s$orig.ident))
      substring <- sub(paste0(unique(s$cellline)),"",substring)
      substring <- sub("ORG1","",substring)
      substring <- sub("ORG2","",substring)
      substring <- sub("-","",substring)
      substring <- sub("-","",substring)
      s$protocol <- unique(substring)
      
      Idents(s) <- s$replicate
      repl.ls.Seurat <- SplitObject(s)
      for(r in repl.ls.Seurat) {
        print(head(colnames(r)))
        saveRDS(r, paste0(OutDirOrig, unique(r$orig.ident),".Rds"))
        saveRDS(r@assays$RNA$counts, paste0(OutDirOrig, unique(r$orig.ident),"_countmatrix.RDS"))
      }
    }
  } else {
    HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(obj), 
                                                 unmapped.as.na = FALSE, 
                                                 map = NULL, 
                                                 species = "human")
    rownames(obj) <- make.unique(HGNC.updated$Suggested.Symbol)
    
    obj <- CreateSeuratObject(obj, 
                              min.cells = 0, 
                              min.features = 0)
    obj$batch <- lane
    
    demultiplex <- read.csv(paste0(OutDirOrig,lane,"/clusters.named.tsv"), sep = "\t")
    demultiplex <- subset(demultiplex, status = "singlet")
    demultiplex <- demultiplex[toupper(demultiplex$assignment.named) %in% toupper(cell.line),]
    
    obj <- subset(obj, cells = demultiplex$barcode)
    obj$cellline <- toupper(demultiplex$assignment.named)
    
    Idents(obj) <- obj$cellline 
    ls.Seurat <- SplitObject(obj)
    
    for(s in ls.Seurat) {
      cellline.ident <- orig.idents[grepl(unique(s$cellline),orig.idents)]
      s$orig.ident <- cellline.ident
      
      if(unique(substr(s$orig.ident, 1,1)) == 5) {
        s$cellline <- "UOFV_1"
      }
      if(unique(substr(s$orig.ident, 1,1)) == 6) {
        s$cellline <- "176_E1"
      }
      if(unique(substr(s$orig.ident, 1,1)) == 7) {
        s$cellline <- "176_E2"
      }
      
      if(grepl("ORG1",unique(s$orig.ident))) {
        s$replicate <- "ORG1"
      } else {
        s$replicate <- "ORG2"
      }

      substring <- substr(s$orig.ident,5,length(s$orig.ident))
      substring <- sub(paste0(unique(s$cellline)),"",substring)
      substring <- sub("ORG1","",substring)
      substring <- sub("ORG2","",substring)
      substring <- sub("-","",substring)
      substring <- sub("-","",substring)
      s$protocol <- unique(substring)
      
      print(head(colnames(s)))
      saveRDS(s, paste0(OutDirOrig, unique(s$orig.ident),".Rds"))
      saveRDS(s@assays$RNA$counts, paste0(OutDirOrig, unique(s$orig.ident),"_countmatrix.RDS"))
    }
  }
}
