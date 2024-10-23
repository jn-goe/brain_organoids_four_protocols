########### experimental conditions colors in vitro ##################
protocol_colors = c("dorsal" = "seagreen",
                    "ventral" = "firebrick",
                    "midbrain" = "royalblue4",
                    "striatum" = "gold2",
                    "pluripotent" = "grey")

protocol_colors_bright = c("dorsal" = "#00ff11ff",
                           "ventral" = "#ff0000ff",
                           "midbrain" = "#003cffff",
                           "striatum" = "#ffd500ff",
                           "pluripotent" = "grey")

cell.line.repl_colors = c("UOFV_1" = "#c21616ff",
                          "XUJA_2" = "#cd80c9ff",
                          "H1" = "#08af8bff",
                          "H9" = "#78bcdeff",
                          "ROZH_5" = "#702ea2ff",
                          "177" = "#e28c2fff",
                          "178"= "#3b47ceff",
                          "176_E1" = "#167711ff",
                          "176_E2" = "#b6dc1fff")

cell.line_colors = cell.line.repl_colors[c("UOFV_1",
                                           "XUJA_2",
                                           "H1",
                                           "H9",
                                           "ROZH_5",
                                           "177",
                                           "178",
                                           "176_E1")]

names(cell.line_colors) = c("UOFV_1",
                            "XUJA_2",
                            "H1",
                            "H9",
                            "ROZH_5",
                            "177",
                            "178",
                            "176")

growth_colors = c("protocol-driven" = "#00A087B2",
                  "cell line-driven" = "#3C5488B2")

condition_colors = c("Residuals" = "grey",
                     "cell.line" = "firebrick2",
                     "protocol" = "forestgreen",
                     "time" = "orange")

########### cell type colors in vitro ##################
dorsal_forebrain_color_palette_fine_num <-
  colorRampPalette(colors = c("#9ddaa6ff", "#56bf66ff", "#005d01ff"))(9)
ventral_forebrain_color_palette_fine_num <-
  colorRampPalette(colors = c("#e7becbff", "#d55780ff","#80024dff"))(12)
striatum_color_palette_fine_num <- c("#d47045ff","#b41b1bff")
midbrain_color_palette_fine_num <- c("#b8bbe7ff", "#4f46b6ff","#a4c9dcff", "#73c7d4ff","#006a8fff")
diencephalon_lin_color_palette_fine_num <- "#a37728ff"
eye_color_palette_fine_num <-
  colorRampPalette(colors = c("#f7e253ff", "#eba600ff"))(6)
oligo_color_palette_fine_num <- c("#b5e82eff", "#9fb309ff")
choroid_plexus_color_palette_fine_num <- "#e71cbfff"
hindbrain_color_palette_fine_num <-
  colorRampPalette(colors = c("#cd97e0ff", "#8e13baff", "#6113baff"))(8)
stromal_color_palette_fine_num <-
  colorRampPalette(colors = c("grey85", "grey25"))(17)

annotation_fine_num_colors <- list("dorsal forebrain_cycl-48" = dorsal_forebrain_color_palette_fine_num[1],
                                   "dorsal forebrain_PC-10" = dorsal_forebrain_color_palette_fine_num[2],
                                   "dorsal forebrain_PC-13" = dorsal_forebrain_color_palette_fine_num[3],
                                   "dorsal forebrain_IP-7" = dorsal_forebrain_color_palette_fine_num[4],
                                   "dorsal forebrain_im-4" = dorsal_forebrain_color_palette_fine_num[5],
                                   "dorsal forebrain_im-20" = dorsal_forebrain_color_palette_fine_num[6],
                                   "dorsal forebrain_UL-6" = dorsal_forebrain_color_palette_fine_num[7],
                                   "dorsal forebrain_DL-22" = dorsal_forebrain_color_palette_fine_num[8],
                                   "dorsal forebrain_DL-25" = dorsal_forebrain_color_palette_fine_num[9],
                                   
                                   "ventral forebrain_cycl-29" = ventral_forebrain_color_palette_fine_num[1],
                                   "ventral forebrain_cycl-36" = ventral_forebrain_color_palette_fine_num[2],
                                   "ventral forebrain_cycl-44" = ventral_forebrain_color_palette_fine_num[3],
                                   "ventral forebrain_PC-17" = ventral_forebrain_color_palette_fine_num[4],
                                   "ventral forebrain_PC-21" = ventral_forebrain_color_palette_fine_num[5],
                                   "ventral forebrain_PC-23" = ventral_forebrain_color_palette_fine_num[6],
                                   "ventral forebrain_IP-15" = ventral_forebrain_color_palette_fine_num[7],
                                   "ventral forebrain_IP-16" = ventral_forebrain_color_palette_fine_num[8],
                                   "ventral forebrain_CGE LGE-11" = ventral_forebrain_color_palette_fine_num[9],
                                   "ventral forebrain_MGE-0" = ventral_forebrain_color_palette_fine_num[10],
                                   "ventral forebrain_MGE-3" = ventral_forebrain_color_palette_fine_num[11],
                                   "ventral forebrain_MGE-32" = ventral_forebrain_color_palette_fine_num[12],
                                   
                                   "striatum-2" = striatum_color_palette_fine_num[1],
                                   "striatum-58" = striatum_color_palette_fine_num[2],
                                   
                                   "midbrain_PC-28" = midbrain_color_palette_fine_num[1],
                                   "midbrain_PC-49" = midbrain_color_palette_fine_num[2],
                                   "midbrain_PC-53" = midbrain_color_palette_fine_num[3],
                                   "midbrain-5" = midbrain_color_palette_fine_num[4],
                                   "midbrain-37" = midbrain_color_palette_fine_num[5],
                                   
                                   "diencephalon lin.-12" = diencephalon_lin_color_palette_fine_num[1],
                                   
                                   "eye_cycl-60" = eye_color_palette_fine_num[1],
                                   "eye_PC-45" = eye_color_palette_fine_num[2],
                                   "eye_ganglion-39" = eye_color_palette_fine_num[3],
                                   "eye_photoreceptors-26" = eye_color_palette_fine_num[4],
                                   "eye_photoreceptors-40" = eye_color_palette_fine_num[5],
                                   "eye_photoreceptors-50" = eye_color_palette_fine_num[6],
                                   
                                   "oligo_PC-19" = oligo_color_palette_fine_num[1],
                                   "oligo-52" = oligo_color_palette_fine_num[2],
                                   
                                   "choroid plexus-31" = choroid_plexus_color_palette_fine_num[1],
                                   
                                   "hindbrain_PC-14" = hindbrain_color_palette_fine_num[1],
                                   "hindbrain_PC-43" = hindbrain_color_palette_fine_num[2],
                                   "hindbrain_PC-54" = hindbrain_color_palette_fine_num[3],
                                   "hindbrain_PC-55" = hindbrain_color_palette_fine_num[4],
                                   "hindbrain_PC-56" = hindbrain_color_palette_fine_num[5],
                                   "hindbrain-1" = hindbrain_color_palette_fine_num[6],
                                   "hindbrain-18" = hindbrain_color_palette_fine_num[7],
                                   "hindbrain-30" = hindbrain_color_palette_fine_num[8],
                                   
                                   "muscle_cycl-33" = stromal_color_palette_fine_num[1],
                                   "muscle_cycl-34" = stromal_color_palette_fine_num[2],
                                   "muscle_cycl-35" = stromal_color_palette_fine_num[3],
                                   "muscle-9" = stromal_color_palette_fine_num[4],
                                   "muscle-27" = stromal_color_palette_fine_num[5],
                                   "muscle-38" = stromal_color_palette_fine_num[6],
                                   "muscle-41" = stromal_color_palette_fine_num[7],
                                   "muscle-46" = stromal_color_palette_fine_num[8],
                                   "muscle-51" = stromal_color_palette_fine_num[9],
                                   
                                   "stromal_fibroblasts-8" = stromal_color_palette_fine_num[10],
                                   "stromal_fibroblasts-24" = stromal_color_palette_fine_num[11],
                                   "stromal_fibroblasts-42" = stromal_color_palette_fine_num[12],
                                   "stromal_epithelial-57" = stromal_color_palette_fine_num[13],
                                   "stromal_epithelial-59" = stromal_color_palette_fine_num[14],
                                   "stromal_epithelial-61" = stromal_color_palette_fine_num[15],
                                   "stromal_immune-47" = stromal_color_palette_fine_num[16],
                                   "stromal_immune-62" = stromal_color_palette_fine_num[17])

dorsal_forebrain_color_palette_fine <-
  colorRampPalette(colors = c("#9ddaa6ff", "#56bf66ff", "#005d01ff"))(6)
ventral_forebrain_color_palette_fine <-
  colorRampPalette(colors = c("#e7becbff", "#d55780ff","#80024dff"))(5)
striatum_color_palette_fine <- "#b41b1bff"
midbrain_color_palette_fine <- c("#b8bbe7ff","#006a8fff")
diencephalon_lin_color_palette_fine <- "#a37728ff"
eye_color_palette_fine <-
  colorRampPalette(colors = c("#f7e253ff", "#eba600ff"))(4)
oligo_color_palette_fine <- c("#b5e82eff", "#9fb309ff")
choroid_plexus_color_palette_fine <- "#e71cbfff"
hindbrain_color_palette_fine <-
  colorRampPalette(colors = c("#cd97e0ff", "#8e13baff", "#6113baff"))(2)
stromal_color_palette_fine <-
  colorRampPalette(colors = c("grey85", "grey25"))(5)

annotation_fine_colors <- list(
  "dorsal forebrain_cycl" = dorsal_forebrain_color_palette_fine[1],
  "dorsal forebrain_PC" = dorsal_forebrain_color_palette_fine[2],
  "dorsal forebrain_IP" = dorsal_forebrain_color_palette_fine[3],
  "dorsal forebrain_im" = dorsal_forebrain_color_palette_fine[4],
  "dorsal forebrain_UL" = dorsal_forebrain_color_palette_fine[5],
  "dorsal forebrain_DL" = dorsal_forebrain_color_palette_fine[6],
  
  "ventral forebrain_cycl" = ventral_forebrain_color_palette_fine[1],
  "ventral forebrain_PC" = ventral_forebrain_color_palette_fine[2],
  "ventral forebrain_IP" = ventral_forebrain_color_palette_fine[3],
  "ventral forebrain_CGE LGE" = ventral_forebrain_color_palette_fine[4],
  "ventral forebrain_MGE" = ventral_forebrain_color_palette_fine[5],
  
  "striatum" = striatum_color_palette_fine[1],
  
  "midbrain_PC" = midbrain_color_palette_fine_num[1],
  "midbrain" = midbrain_color_palette_fine_num[2],

  "diencephalon lin." = diencephalon_lin_color_palette_fine[1],
  
  "eye_cycl" = eye_color_palette_fine[1],
  "eye_PC" = eye_color_palette_fine[2],
  "eye_ganglion" = eye_color_palette_fine[3],
  "eye_photoreceptors" = eye_color_palette_fine[4],
  
  "oligo_PC" = oligo_color_palette_fine[1],
  "oligo" = oligo_color_palette_fine[2],
  
  "choroid plexus" = choroid_plexus_color_palette_fine[1],
  
  "hindbrain_PC" = hindbrain_color_palette_fine[1],
  "hindbrain" = hindbrain_color_palette_fine[2],
  
  "muscle_cycl" = stromal_color_palette_fine[1],
  "muscle" = stromal_color_palette_fine[2],
  "stromal_fibroblasts" = stromal_color_palette_fine[3],
  "stromal_epithelial" = stromal_color_palette_fine[4],
  "stromal_immune" = stromal_color_palette_fine[5]
)

dorsal_forebrain_color_palette_coarse <- c("#56bf66ff", "#005d01ff")
ventral_forebrain_color_palette_coarse <- c("#d55780ff","#80024dff")
striatum_color_palette_coarse<- "#b41b1bff"
midbrain_color_palette_coarse <- c("#b8bbe7ff","#006a8fff")
diencephalon_lin_color_palette_coarse <- "#a37728ff"
eye_color_palette_coarse <- "#e9c736ff"
oligo_color_palette_coarse <- "#b5e82eff"
choroid_plexus_color_palette_coarse <- "#e71cbfff"
hindbrain_color_palette_coarse <- c("#cd97e0ff", "#8e13baff")
stromal_color_palette_coarse <- c("grey35","grey75")

annotation_coarse_colors <- list(
  "dorsal forebrain prog." = dorsal_forebrain_color_palette_coarse[1],
  "dorsal forebrain" = dorsal_forebrain_color_palette_coarse[2],
  
  "ventral forebrain prog." = ventral_forebrain_color_palette_coarse[1],
  "ventral forebrain" = ventral_forebrain_color_palette_coarse[2],
  
  "striatum" = striatum_color_palette_coarse[1],
  
  "midbrain prog." = midbrain_color_palette_coarse[1],
  "midbrain" = midbrain_color_palette_coarse[2],
  
  "diencephalon lin." = diencephalon_lin_color_palette_coarse[1],
  
  "eye" = eye_color_palette_coarse[1],
  
  "oligo" = oligo_color_palette_coarse[1],
  
  "choroid plexus" = choroid_plexus_color_palette_coarse[1],
  
  "hindbrain prog." = hindbrain_color_palette_coarse[1],
  "hindbrain" = hindbrain_color_palette_coarse[2],
  
  "muscle" = stromal_color_palette_coarse[1],
  "stromal"= stromal_color_palette_coarse[2])

########### cell type colors in vivo braun ##################

Region_colors <-list("telencephalon" = "#56bf66ff",
                     "midbrain" = "#73c7d4ff",
                     "diencephalon" = "#d55780ff",
                     "cerebellum" = "#8e13baff")

Subregion_colors <- list("cortex" = "#56bf66ff",
                         "hippocampus" = "#c0ddb4ff",     
                         "striatum" = "#b41b1bff",         
                         "dorsal midbrain" = "#73c7d4ff",  
                         "ventral midbrain"= "#b8bbe7ff",
                         "thalamus" = "#bfb22fff",
                         "hypothalamus" = "#a37728ff",     
                         "cerebellum" = "#8e13baff")

CellClass_colors <- list("neural crest" = "#d2e2d9fe",
                         "glioblast" = "#a8dfd0fe", 
                         "radial glia" = "#82d1c9fe",
                         "neuroblast" = "#3dbea8fe",
                         "neuronal ipc" = "#2cc7cffe", 
                         "neuron" = "#3892a8fe",
                         "oligo" = "#9fb309ff",  
                         "fibroblast" = "grey60",
                         "immune" = "grey80",
                         "vascular" = "grey20",
                         "erythrocyte" = "grey40")

cortex_color_palette_braun <- colorRampPalette(colors = c("#c0ddb4ff", "#56bf66ff","#005d01ff"))(5)
hippocampus_color_palette_braun <- colorRampPalette(colors = c("#f1daa8ff","#d7af03ff"))(5)
striatum_color_palette_braun <- colorRampPalette(colors = c("#d47045ff","#b41b1bff"))(5)
midbrain_dorsal_color_palette_braun <- colorRampPalette(colors = c("#a4c9dcff", "#73c7d4ff","#006a8fff"))(5)
midbrain_ventral_color_palette_braun <- colorRampPalette(colors = c("#b8bbe7ff", "#4f46b6ff"))(5)
thalamus_color_palette_braun <- colorRampPalette(colors = c("#dae29eff", "#c5d733ff","#bfb22fff","#a37728ff"))(10)[1:5]
hypothalamus_color_palette_braun <- colorRampPalette(colors = c("#dae29eff", "#c5d733ff","#bfb22fff","#a37728ff"))(10)[6:10]
cerebellum_color_palette_braun <- colorRampPalette(colors = c("#cd97e0ff", "#8e13baff", "#6113baff"))(5)

annotation_braun_colors <- list("neural crest" = "#d2e2d9fe",
                          
                          "cortex_glioblast" = cortex_color_palette_braun[1],
                          "cortex_radial glia" = cortex_color_palette_braun[2],
                          "cortex_neuroblast" = cortex_color_palette_braun[3],
                          "cortex_neuronal ipc" = cortex_color_palette_braun[4],
                          "cortex_neuron" = cortex_color_palette_braun[5],
                          
                          "hippocampus_glioblast" = hippocampus_color_palette_braun[1],
                          "hippocampus_radial glia" = hippocampus_color_palette_braun[2],
                          "hippocampus_neuroblast" = hippocampus_color_palette_braun[3],
                          "hippocampus_neuronal ipc" = hippocampus_color_palette_braun[4],
                          "hippocampus_neuron" = hippocampus_color_palette_braun[5],
                          
                          "striatum_glioblast" = striatum_color_palette_braun[1],
                          "striatum_radial glia" = striatum_color_palette_braun[2],
                          "striatum_neuroblast" = striatum_color_palette_braun[3],
                          "striatum_neuronal ipc" = striatum_color_palette_braun[4],
                          "striatum_neuron" = striatum_color_palette_braun[5],
                          
                          "dorsal midbrain_glioblast" = midbrain_dorsal_color_palette_braun[1],
                          "dorsal midbrain_radial glia" = midbrain_dorsal_color_palette_braun[2],
                          "dorsal midbrain_neuroblast" = midbrain_dorsal_color_palette_braun[3],
                          "dorsal midbrain_neuronal ipc" = midbrain_dorsal_color_palette_braun[4],
                          "dorsal midbrain_neuron" = midbrain_dorsal_color_palette_braun[5],
                          
                          "ventral midbrain_glioblast" = midbrain_ventral_color_palette_braun[1],
                          "ventral midbrain_radial glia" = midbrain_ventral_color_palette_braun[2],
                          "ventral midbrain_neuroblast" = midbrain_ventral_color_palette_braun[3],
                          "ventral midbrain_neuronal ipc" = midbrain_ventral_color_palette_braun[4],
                          "ventral midbrain_neuron" = midbrain_ventral_color_palette_braun[5],
                          
                          "thalamus_glioblast" = thalamus_color_palette_braun[1],
                          "thalamus_radial glia" = thalamus_color_palette_braun[2],
                          "thalamus_neuroblast" = thalamus_color_palette_braun[3],
                          "thalamus_neuronal ipc" = thalamus_color_palette_braun[4],
                          "thalamus_neuron" = thalamus_color_palette_braun[5],
                          
                          "hypothalamus_glioblast" = hypothalamus_color_palette_braun[1],
                          "hypothalamus_radial glia" = hypothalamus_color_palette_braun[2],
                          "hypothalamus_neuroblast" = hypothalamus_color_palette_braun[3],
                          "hypothalamus_neuronal ipc" = hypothalamus_color_palette_braun[4],
                          "hypothalamus_neuron" = hypothalamus_color_palette_braun[5],
                          
                          "cerebellum_glioblast" = cerebellum_color_palette_braun[1],
                          "cerebellum_radial glia" = cerebellum_color_palette_braun[2],
                          "cerebellum_neuroblast" = cerebellum_color_palette_braun[3],
                          "cerebellum_neuronal ipc" = cerebellum_color_palette_braun[4],
                          "cerebellum_neuron" = cerebellum_color_palette_braun[5],
                          
                          "oligo" = "#9fb309ff",  
                          "fibroblast" = "grey60",
                          "immune" = "grey80",
                          "vascular" = "grey20",
                          "erythrocyte" = "grey40")

########### cell type colors in vivo bhaduri ##################

neocortex_color_palette_bhaduri <- colorRampPalette(colors = c("#c0ddb4ff", "#56bf66ff","#005d01ff"))(10)[1:5]
allocortex_color_palette_bhaduri <- colorRampPalette(colors = c("#c0ddb4ff", "#56bf66ff","#005d01ff"))(10)[6:10]
ge_color_palette_bhaduri <- colorRampPalette(colors = c("#e7becbff", "#d55780ff","#80024dff"))(9)[1:5]
claustrum_color_palette_bhaduri <- colorRampPalette(colors = c("#e7becbff", "#d55780ff","#80024dff"))(9)[6:9]
striatum_color_palette_bhaduri <- colorRampPalette(colors = c("#d47045ff","#b41b1bff"))(4)
thalamus_color_palette_bhaduri <- colorRampPalette(colors = c("#dae29eff", "#c5d733ff","#bfb22fff","#a37728ff"))(8)[1:4]
hypothalamus_color_palette_bhaduri <- colorRampPalette(colors = c("#dae29eff", "#c5d733ff","#bfb22fff","#a37728ff"))(8)[5:8]

annotation_bhaduri_colors <- list("neocortex_dividing" = neocortex_color_palette_bhaduri[1],
                          "neocortex_rg" = neocortex_color_palette_bhaduri[2],
                          "neocortex_ipc" = neocortex_color_palette_bhaduri[3],
                          "neocortex_neuron" = neocortex_color_palette_bhaduri[4],
                          "neocortex_interneuron" = neocortex_color_palette_bhaduri[5],
                          
                          "allocortex_dividing" = allocortex_color_palette_bhaduri[1],
                          "allocortex_rg" = allocortex_color_palette_bhaduri[2],
                          "allocortex_ipc" = allocortex_color_palette_bhaduri[3],
                          "allocortex_neuron" = allocortex_color_palette_bhaduri[4],
                          "allocortex_interneuron" = allocortex_color_palette_bhaduri[5],
                          
                          "claustrum_dividing" = claustrum_color_palette_bhaduri[1],
                          "claustrum_rg" = claustrum_color_palette_bhaduri[2],
                          "claustrum_interneuron" = claustrum_color_palette_bhaduri[3],
                          "claustrum_neuron" = claustrum_color_palette_bhaduri[4],
                          
                          "ge_dividing" = ge_color_palette_bhaduri[1],
                          "ge_rg" = ge_color_palette_bhaduri[2],
                          "ge_ipc" = ge_color_palette_bhaduri[3],
                          "ge_neuron" = ge_color_palette_bhaduri[4],
                          "ge_interneuron" = ge_color_palette_bhaduri[5],
                          
                          "hypothalamus_dividing" = hypothalamus_color_palette_bhaduri[1],
                          "hypothalamus_rg" = hypothalamus_color_palette_bhaduri[2],
                          "hypothalamus_ipc" = hypothalamus_color_palette_bhaduri[3],
                          "hypothalamus_neuron" = hypothalamus_color_palette_bhaduri[4],
                          
                          "thalamus_dividing" = thalamus_color_palette_bhaduri[1],
                          "thalamus_rg" = thalamus_color_palette_bhaduri[2],
                          "thalamus_interneuron" = thalamus_color_palette_bhaduri[3],
                          "thalamus_neuron" = thalamus_color_palette_bhaduri[4],
                          
                          "striatum_dividing" = striatum_color_palette_bhaduri[1],
                          "striatum_rg" = striatum_color_palette_bhaduri[2],
                          "striatum_interneuron" = striatum_color_palette_bhaduri[3],
                          "striatum_neuron" = striatum_color_palette_bhaduri[4],
                          
                          "oligo" = "#add768ff",
                          "astrocyte" = "#d0dbb0ff",
                          
                          "microglia" = "grey90",
                          
                          "endo" = "grey40",
                          "outlier" = "black",
                          "vascular" = "grey70")