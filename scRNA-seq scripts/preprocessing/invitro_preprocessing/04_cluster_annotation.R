rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
library(Seurat)
library(ggplot2)

source("./helper_functions/colors_metadata.R")

# Setup ----------------------------
OutDirOrig <- "./scRNAseq_analysis/allinone/"
setwd(OutDirOrig)

DataDir <- "./scRNAseq_analysis/data/"
combined.obj <- readRDS(paste0(DataDir,"combined.obj.RDS"))

combined.obj <- FindClusters(combined.obj, resolution = 3.5)

######################### ANNOTATION #########################

annotation_fine_num <- list("0" = "ventral forebrain_MGE-0",
                            "1" = "hindbrain-1",
                            "2" = "striatum-2",
                            "3" = "ventral forebrain_MGE-3",
                            "4" = "dorsal forebrain_im-4",
                            "5" = "midbrain-5",
                            "6" = "dorsal forebrain_UL-6",
                            "7" = "dorsal forebrain_IP-7",
                            "8" = "stromal_fibroblasts-8",
                            "9" = "muscle-9",
                            "10" = "dorsal forebrain_PC-10",
                            "11" = "ventral forebrain_CGE LGE-11",
                            "12" = "diencephalon lin.-12",
                            "13" = "dorsal forebrain_PC-13",
                            "14" = "hindbrain_PC-14",
                            "15" = "ventral forebrain_IP-15",
                            "16" = "ventral forebrain_IP-16",
                            "17" = "ventral forebrain_PC-17",
                            "18" = "hindbrain-18",
                            "19" = "oligo_PC-19",
                            "20" = "dorsal forebrain_im-20",
                            "21" = "ventral forebrain_PC-21",
                            "22" = "dorsal forebrain_DL-22",
                            "23" = "ventral forebrain_PC-23",
                            "24" = "stromal_fibroblasts-24",
                            "25" = "dorsal forebrain_DL-25",
                            "26" = "eye_photoreceptors-26",
                            "27" = "muscle-27",
                            "28" = "midbrain_PC-28", 
                            "29" = "ventral forebrain_cycl-29",
                            "30" = "hindbrain-30",
                            "31" = "choroid plexus-31",
                            "32" = "ventral forebrain_MGE-32",
                            "33" = "muscle_cycl-33",
                            "34" = "muscle_cycl-34",
                            "35" = "muscle_cycl-35",
                            "36" = "ventral forebrain_cycl-36",
                            "37" = "midbrain-37", 
                            "38" = "muscle-38",
                            "39" = "eye_ganglion-39",
                            "40" = "eye_photoreceptors-40",
                            "41" = "muscle-41",
                            "42" = "stromal_fibroblasts-42",
                            "43" = "hindbrain_PC-43", 
                            "44" = "ventral forebrain_cycl-44",
                            "45" = "eye_PC-45",
                            "46" = "muscle-46",
                            "47" = "stromal_immune-47", 
                            "48" = "dorsal forebrain_cycl-48",
                            "49" = "midbrain_PC-49",
                            "50" = "eye_photoreceptors-50",
                            "51" = "muscle-51",
                            "52" = "oligo-52",
                            "53" = "midbrain_PC-53",
                            "54" = "hindbrain_PC-54",
                            "55" = "hindbrain_PC-55",
                            "56" = "hindbrain_PC-56",
                            "57" = "stromal_epithelial-57",
                            "58" = "striatum-58",
                            "59" = "stromal_epithelial-59",
                            "60" = "eye_cycl-60",
                            "61" = "stromal_epithelial-61", 
                            "62" = "stromal_immune-62") 

annotation_fine <- unlist(lapply(annotation_fine_num, function(x) strsplit(x,"-")[[1]][1]))
annotation_coarse <- annotation_fine

annotation_coarse[grepl("_cycl|_PC|_IP", annotation_coarse)] <- 
  sub("_cycl|_PC|_IP", " prog.", 
      annotation_coarse[grepl("_cycl|_PC|_IP", annotation_coarse)])
annotation_coarse[grepl("_MGE|_CGE LGE|_UL|_DL|_im|_photoreceptors|_fibroblasts|_ganglion|_unspecified", annotation_coarse)] <- 
  sub("_MGE|_CGE LGE|_UL|_DL|_im|_photoreceptors|_fibroblasts|_ganglion|_unspecified", "", 
      annotation_coarse[grepl("_MGE|_CGE LGE|_UL|_DL|_im|_photoreceptors|_fibroblasts|_ganglion|_unspecified", annotation_coarse)])
annotation_coarse[grepl("stromal", annotation_coarse)] <- "stromal"
annotation_coarse[grepl("muscle", annotation_coarse)] <- "muscle"
annotation_coarse[grepl("eye", annotation_coarse)] <- "eye"
annotation_coarse[grepl("oligo", annotation_coarse)] <- "oligo"
annotation_coarse[grepl("midbrain prog.", annotation_coarse)] <- "midbrain prog."

annotation_fine_num <- factor(annotation_fine_num, levels = 
                                c(
                                  "dorsal forebrain_cycl-48",
                                  "dorsal forebrain_PC-10",
                                  "dorsal forebrain_PC-13",
                                  "dorsal forebrain_IP-7",
                                  "dorsal forebrain_im-4",
                                  "dorsal forebrain_im-20",
                                  "dorsal forebrain_UL-6",
                                  "dorsal forebrain_DL-22",
                                  "dorsal forebrain_DL-25",
                                  
                                  "ventral forebrain_cycl-29",
                                  "ventral forebrain_cycl-36", 
                                  "ventral forebrain_cycl-44",
                                  "ventral forebrain_PC-17",
                                  "ventral forebrain_PC-21",
                                  "ventral forebrain_PC-23",
                                  "ventral forebrain_IP-15",
                                  "ventral forebrain_IP-16",
                                  "ventral forebrain_CGE LGE-11",
                                  "ventral forebrain_MGE-0",
                                  "ventral forebrain_MGE-3",
                                  "ventral forebrain_MGE-32",

                                  "striatum-2",
                                  "striatum-58",
                                  
                                  "midbrain_PC-28",
                                  "midbrain_PC-49",
                                  "midbrain_PC-53",
                                  "midbrain-5",
                                  "midbrain-37",
                                  
                                  "diencephalon lin.-12",

                                  "eye_cycl-60",
                                  "eye_PC-45",
                                  "eye_ganglion-39",
                                  "eye_photoreceptors-26",
                                  "eye_photoreceptors-40",
                                  "eye_photoreceptors-50",

                                  "oligo_PC-19",
                                  "oligo-52",

                                  "choroid plexus-31",
                                  
                                  "hindbrain_PC-14",
                                  "hindbrain_PC-43",
                                  "hindbrain_PC-54",
                                  "hindbrain_PC-55",
                                  "hindbrain_PC-56",
                                  "hindbrain-1",
                                  "hindbrain-18",
                                  "hindbrain-30",

                                  "muscle_cycl-33",
                                  "muscle_cycl-34",
                                  "muscle_cycl-35",
                                  "muscle-9",
                                  "muscle-27",
                                  "muscle-38",
                                  "muscle-41",
                                  "muscle-46",
                                  "muscle-51",

                                  "stromal_fibroblasts-8",
                                  "stromal_fibroblasts-24",
                                  "stromal_fibroblasts-42",
                                  "stromal_epithelial-57",
                                  "stromal_epithelial-59",
                                  "stromal_epithelial-61",
                                  "stromal_immune-47",
                                  "stromal_immune-62"))

combined.obj@misc$annotation_fine_num_colors <- annotation_fine_num_colors

Idents(combined.obj) <- combined.obj$RNA_snn_res.3.5
combined.obj <- RenameIdents(combined.obj, annotation_fine_num)
combined.obj$annotation_fine_num <- Idents(combined.obj) 

annotation_fine <- factor(annotation_fine, levels =
                                c(
                                  "dorsal forebrain_cycl",
                                  "dorsal forebrain_PC",
                                  "dorsal forebrain_IP",
                                  "dorsal forebrain_im",
                                  "dorsal forebrain_UL",
                                  "dorsal forebrain_DL",

                                  "ventral forebrain_cycl",
                                  "ventral forebrain_PC",
                                  "ventral forebrain_IP",
                                  "ventral forebrain_CGE LGE",
                                  "ventral forebrain_MGE",

                                  "striatum",

                                  "midbrain_PC",
                                  "midbrain",

                                  "diencephalon lin.",

                                  "eye_cycl",
                                  "eye_PC",
                                  "eye_ganglion",
                                  "eye_photoreceptors",

                                  "oligo_PC",
                                  "oligo",

                                  "choroid plexus",

                                  "hindbrain_PC",
                                  "hindbrain",

                                  "muscle_cycl",
                                  "muscle",
                                  "stromal_fibroblasts",
                                  "stromal_epithelial",
                                  "stromal_immune"
                                ))
 
combined.obj@misc$annotation_fine_colors <- annotation_fine_colors

Idents(combined.obj) <- combined.obj$RNA_snn_res.3.5
combined.obj <- RenameIdents(combined.obj, annotation_fine)
combined.obj$annotation_fine <- Idents(combined.obj) 

annotation_coarse <- factor(annotation_coarse, levels =
                            c("dorsal forebrain prog.",
                              "dorsal forebrain" ,

                              "ventral forebrain prog.",
                              "ventral forebrain",

                              "striatum",

                              "midbrain prog.",
                              "midbrain",
                              
                              "diencephalon lin.",

                              "eye",

                              "oligo",

                              "choroid plexus",

                              "hindbrain prog.",
                              "hindbrain",

                              "muscle",
                              "stromal"
           ))

combined.obj@misc$annotation_coarse_colors <- annotation_coarse_colors

Idents(combined.obj) <- combined.obj$RNA_snn_res.3.5
combined.obj <- RenameIdents(combined.obj, annotation_coarse)
combined.obj$annotation_coarse <- Idents(combined.obj)

saveRDS(combined.obj, paste0(DataDir,"combined.obj.RDS"))
