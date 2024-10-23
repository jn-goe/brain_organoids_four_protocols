plot_QC <- function(s_obj
                    , nfeature.range = c(250,max(s_obj$nFeature_RNA))
                    , ncount.range = c(500,max(s_obj$nCount_RNA))
                    , mito.range = c(0.001,0.2)
                    , ribo.range = c(0.001,0.5))  {

    mm <- s_obj@meta.data

    below.nFeature_RNA <- nfeature.range[1]
    above.nFeature_RNA <- nfeature.range[2]

    below.nCount_RNA <- ncount.range[1]
    above.nCount_RNA <- ncount.range[2]

    below.mito <- mito.range[1]
    above.mito <- mito.range[2]

    below.ribo <- ribo.range[1]
    above.ribo <- ribo.range[2]

    mm$colour.thr.nFeature <- cut(mm$nFeature_RNA,
                                  breaks = c(-Inf, above.nFeature_RNA, below.nFeature_RNA, Inf),
                                  labels = c("coral","skyblue","coral"))

    mm$colour.thr.nCount_RNA <- cut(mm$nCount_RNA,
                                    breaks = c(-Inf, above.nCount_RNA, below.nCount_RNA, Inf),
                                    labels = c("coral","skyblue","coral"))

    mm$colour.thr.mito <- cut(mm$percent.mito,
                              breaks = c(-Inf, above.mito, below.mito, Inf),
                              labels = c("coral","skyblue","coral"))

    mm$colour.thr.ribo <- cut(mm$percent.ribo,
                              breaks = c(-Inf, above.ribo, below.ribo, Inf),
                              labels = c("coral","skyblue","coral"))

    filt.nFeature_RNA <- (mm$nFeature_RNA > below.nFeature_RNA & mm$nFeature_RNA < above.nFeature_RNA)
    filt.nCount_RNA <- (mm$nCount_RNA > below.nCount_RNA & mm$nCount_RNA < above.nCount_RNA)
    filt.mito <- (mm$percent.mito > below.mito & mm$percent.mito < above.mito)
    filt.ribo <- (mm$percent.ribo > below.ribo & mm$percent.ribo < above.ribo)

    mm <- data.frame(nFeature_RNA = mm$nFeature_RNA,
                     nCount_RNA = mm$nCount_RNA,
                     percent.mito = mm$percent.mito,
                     percent.ribo = mm$percent.ribo,
                     filt.nFeature_RNA, filt.mito,
                     filt.ribo, filt.nCount_RNA,
                     colour.thr.nFeature = mm$colour.thr.nFeature,
                     colour.thr.nCount_RNA = mm$colour.thr.nCount_RNA,
                     colour.thr.mito = mm$colour.thr.mito,
                     colour.thr.ribo = mm$colour.thr.ribo
    )

    A <- ggplot2::ggplot(mm, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA)) +
    ggplot2::geom_point(alpha = 0.25, size =  0.75,  show.legend = FALSE,
                        ggplot2::aes(color = filt.nFeature_RNA & filt.nCount_RNA & filt.mito & filt.ribo)  ) +
    ggplot2::scale_y_log10() + ggplot2::scale_x_log10() +
    ggplot2::annotation_logticks()

    B <- ggplot2::ggplot(mm, ggplot2::aes(x = percent.ribo, y = percent.mito)) +
    ggplot2::geom_point(alpha = 0.25, size =  0.75,  show.legend = FALSE,
                        ggplot2::aes(color = filt.nFeature_RNA & filt.nCount_RNA & filt.mito & filt.ribo)) +
    ggplot2::scale_y_log10() + #ggplot2::scale_x_log10() +
    ggplot2::annotation_logticks() +
    ggplot2::geom_hline(yintercept = below.mito) +
    ggplot2::geom_hline(yintercept = above.mito) +
    ggplot2::geom_vline(xintercept = below.ribo) +
    ggplot2::geom_vline(xintercept = above.ribo)

    C <- ggplot2::ggplot(data = mm, ggplot2::aes(x = nFeature_RNA, fill = colour.thr.nFeature)) +
    ggplot2::geom_histogram(binwidth = (2*IQR((mm$nFeature_RNA))/length(mm$nFeature_RNA)^(1/3))) +
    Seurat::NoLegend() +
    ggplot2::geom_vline(xintercept = below.nFeature_RNA) +
    ggplot2::geom_vline(xintercept = above.nFeature_RNA)

    D <- ggplot2::ggplot(data = mm, ggplot2::aes(x = nCount_RNA, fill = colour.thr.nCount_RNA)) +
    ggplot2::geom_histogram(binwidth = (2*IQR((mm$nCount_RNA))/length(mm$nCount_RNA)^(1/3))) +
    Seurat::NoLegend() +
    ggplot2::geom_vline(xintercept = below.nCount_RNA) +
    ggplot2::geom_vline(xintercept = above.nCount_RNA)

    E <- ggplot2::ggplot(data = mm, ggplot2::aes(x = percent.mito, fill = colour.thr.mito)) +
    ggplot2::geom_histogram(binwidth = (2*IQR(log(mm$percent.mito))/length(mm$percent.mito)^(1/3))) +
    Seurat::NoLegend() +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = below.mito) +
    ggplot2::geom_vline(xintercept = above.mito)

    G <- ggplot2::ggplot(data = mm, ggplot2::aes(x = percent.ribo, fill = colour.thr.ribo)) +
    ggplot2::geom_histogram(binwidth = (2*IQR(mm$percent.ribo)/length(mm$percent.ribo)^(1/3))) +
    Seurat::NoLegend() +
    ggplot2::geom_vline(xintercept = below.ribo) +
    ggplot2::geom_vline(xintercept = above.ribo)

    H <- ggplot2::ggplot(mm, ggplot2::aes(x = nFeature_RNA, y = log(percent.mito))) +
    ggplot2::geom_point(alpha = 0.25, size = 0.75, show.legend = FALSE,
                        ggplot2::aes(color = filt.nFeature_RNA & filt.nCount_RNA & filt.mito & filt.ribo)) +
    ggplot2::geom_hline(yintercept = log(below.mito)) +
    ggplot2::geom_hline(yintercept = log(above.mito)) +
    ggplot2::geom_vline(xintercept = below.nFeature_RNA) +
    ggplot2::geom_vline(xintercept = above.nFeature_RNA)

    I <- ggplot2::ggplot(mm, ggplot2::aes(x = nFeature_RNA, y = percent.ribo)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.75, show.legend = FALSE,
                        ggplot2::aes(color = filt.nFeature_RNA & filt.nCount_RNA & filt.mito & filt.ribo)) +
    ggplot2::annotation_logticks() +
    ggplot2::geom_hline(yintercept = below.ribo) +
    ggplot2::geom_hline(yintercept = above.ribo) +
    ggplot2::geom_vline(xintercept = below.nFeature_RNA) +
    ggplot2::geom_vline(xintercept = above.nFeature_RNA)

    J1 <- Seurat::VlnPlot(s_obj, features = c("nFeature_RNA"))
    J2 <- Seurat::VlnPlot(s_obj, features = c("nCount_RNA"))
    J3 <- Seurat::VlnPlot(s_obj, features = c("percent.mito"))

    plotlist <- list(C,D,A,
                     B,E,G)

    s_obj$passQC <- filt.nFeature_RNA & filt.nCount_RNA & filt.mito & filt.ribo

    return(list(plotlist,s_obj))
}

automatic_QC <- function(obj, 
                         mito.range = c(0.001,0.15), 
                         ribo.range = c(0.001,0.5), 
                         nCount_RNA.range = c(500, 20000),
                         nFeature_RNA.range = c(500, 20000), 
                         normal_based = F) {
  
  Cell.QC.Stat <- obj@meta.data
  
  if(!normal_based) { 
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.mito < mito.range[2]) %>% filter(percent.mito > mito.range[1])
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.ribo < ribo.range[2]) %>% filter(percent.ribo > ribo.range[1])
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(nCount_RNA < nCount_RNA.range[2]) %>% filter(nCount_RNA > nCount_RNA.range[1])
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(nFeature_RNA < nFeature_RNA.range[2]) %>% filter(nFeature_RNA > nFeature_RNA.range[1])
  }
  
  if(normal_based) {
    Cell.QC.Stat$percent.mito <- abs(scale((Cell.QC.Stat$percent.mito)))
    Cell.QC.Stat$percent.ribo <- abs(scale((Cell.QC.Stat$percent.ribo)))
    Cell.QC.Stat$nCount_RNA <- abs(scale((Cell.QC.Stat$nCount_RNA)))
    Cell.QC.Stat$nFeature_RNA <- abs(scale((Cell.QC.Stat$nFeature_RNA)))
    
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.mito < qnorm(0.975)) 
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.ribo < qnorm(0.975)) 
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(nCount_RNA < qnorm(0.975)) 
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(nFeature_RNA < qnorm(0.975)) 
  }
  
  keep.cells <- rownames(Cell.QC.Stat)
  return(keep.cells)
}

fun_max_ratio <- function(x) {
  res <- max(x)/sum(x)
  if(is.nan(res)) {
    return(0)
  } else {
    return(res)
  }
}

fun_max_which <- function(x) {
  return(names(which.max(x)))
}
