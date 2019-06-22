
## Dependencies

library(GenomicRanges)
library(bumphunter)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(RColorBrewer)
library(derfinder)
library(derfinderPlot)
library(gridExtra)
library(cowplot)


## Function 4: Make four plots of nuclear, NeuN-sorted, degradation, and cell type-specific expression of REGION

four_panels = function(REGION, PATH="Default", JUNCTIONS="FALSE") {
  
  load("./data/probe_design_four_plot_phenotype_data.rda")
  load("./data/region_info_objects.rda")
  
  gr = GRanges(REGION)
  nearestAnnotation = matchGenes(x = gr, subject = genes)
  
  pd = list(Sep = pdSep, Deg = pdDeg, Cell = pdCell, Sort = pdSort)
  
  files = list(Sep = rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Sep', samplepatt='*.bw$', fileterm=NULL),
               Deg = rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Deg', samplepatt='*.bw$', fileterm=NULL),
               Cell = rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Cell', samplepatt='*.bw$', fileterm=NULL),
               Sort = rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Sort', samplepatt='*.bw$', fileterm=NULL))
  files = mapply(function(bw, mpd) bw[match(mpd$BigWig, names(bw))], files, pd, SIMPLIFY = F)
  
  regionCov = mapply(function(bw, mpd) getRegionCoverage(regions=gr, totalMapped = mpd$sumMapped,
                                                         files = bw), files, pd, SIMPLIFY = F)
  
  if (JUNCTIONS!="FALSE") {
    
    regionCov = lapply(regionCov, function(x) list(do.call(rbind, x)))
    covMat = lapply(regionCov, function(x) t(sapply(x, colSums)/100))
    
    bg = lapply(covMat, function(x) matrix(sum(width(gr)), nc = ncol(x), nr = nrow(x))/1e3)
    
    coords = paste0(seqnames(gr)[1], ":", min(start(gr)),"-", max(end(gr)), ":", strand(gr)[1])
    mains = paste0(coords, " (", sum(width(gr)), " bp)\n", nearestAnnotation$name)
    
  } else {
    
    covMat = lapply(regionCov, function(x) t(sapply(x, colSums)/100))
    
    bg = lapply(covMat, function(x) matrix(width(gr), nc = ncol(x), nr = nrow(x))/1e3)

    coords = paste0(seqnames(gr), ":", start(gr),"-", end(gr), ":", strand(gr))
    mains = paste0(coords, " (", width(gr), " bp)\n", nearestAnnotation$name)

  }
  
  covMat = mapply(function(cv, b) cv/b, covMat, bg, SIMPLIFY = F)
  covMat = lapply(covMat, function(x) as.matrix(log2(x + 1)))
  theRanges = range(do.call(cbind, covMat))
  
  if (PATH=="Default") {
    
    pdf("./four_panels.pdf", h = 12, w = 9)
  
  } else {
    
    pdf(paste0(PATH, "four_panels.pdf"), h = 12, w = 9)
    
  }
  
  for(j in 1:length(regionCov[[1]])) {
    
    SepDat = data.frame(Cov = covMat$Sep[j,], LabelFrac = pdSep$LabelFrac, Shortlabels = pdSep$Shortlabels)
    Sep = ggplot(SepDat, aes(x = LabelFrac, y = Cov)) +
      theme_bw() +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 2, aes(fill = Shortlabels), pch=21, color="black") +
      scale_fill_manual(values = brewer.pal(8,"Paired")) +
      labs(fill="") + 
      ylim(theRanges) +
      ylab("Log2(Adj Read/kb)") + xlab("") +
      theme(text = element_text(size = 20),
            legend.text = element_text(size = 12),
            legend.position=c(0.5, 0.1),
            legend.background = element_rect(fill = "transparent", linetype="solid", colour ="black"), 
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      guides(fill=guide_legend(ncol=4)) +
      ggtitle("Separation")
    
    
    DegDat = data.frame(Cov = covMat$Deg[j,], DegradationTime = pdDeg$DegradationTime, 
                        LibraryProtocol = pdDeg$LibraryProtocol, BrNum = pdDeg$BrNum)
    Deg = ggplot(DegDat, aes(x = DegradationTime, y = Cov, fill = BrNum, color = BrNum)) +
      theme_bw() +
      geom_line(aes(linetype = LibraryProtocol), lwd = 1.3) +
      geom_point(cex = 2, pch = 21, aes(fill = BrNum), color="black") +
      scale_color_manual(values = brewer.pal(5,"Dark2")) +
      scale_fill_manual(values = brewer.pal(5,"Dark2")) +
      ylim(theRanges) +
      ylab("") + xlab("") +
      theme(text = element_text(size = 20),
            legend.text = element_text(size = 12),
            legend.position=c(0.5, 0.1),
            legend.background = element_rect(fill = "transparent", linetype="solid", colour ="black"), 
            legend.title = element_blank()) +
      guides(linetype=guide_legend(ncol=2), fill = FALSE, color = FALSE) +
      ggtitle("Degradation")
    
    
    SortDat = data.frame(Cov = covMat$Sort[j,], Label = gsub(":", "\n", levels(factor(pdSort$Label))))
    Sort = ggplot(SortDat, aes(x = Label, y = Cov)) +
      theme_bw() +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 2, aes(fill = Label), pch=21, color="black") +
      scale_fill_manual(values = unique(pdSort$col)) +
      ylim(theRanges) +
      ylab("Log2(Adj Read/kb)") + xlab("") +
      theme(text = element_text(size = 20),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      guides(fill=FALSE) +
      ggtitle("Sorted")
    
    
    CellDat = data.frame(Cov = covMat$Cell[j,], Cell_type = pdCell$Cell_type)
    Cell = ggplot(CellDat, aes(x = Cell_type, y = Cov)) +
      theme_bw() +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(fill = Cell_type), pch=21, color="black") +
      scale_fill_manual(values = c("#FFBB78","#FF9896","#C5B0D5","#9467BD","#1F77B4",
                                   "#2CA02C","#D62728","#FF7F0E","#AEC7E8","#98DF8A")) +
      ylim(theRanges) +
      ylab("") + xlab("") +
      theme(text = element_text(size = 20),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      guides(fill=FALSE) +
      ggtitle("Single Cells")
    
    p = plot_grid(Sep, Deg, Sort, Cell, ncol = 2, align = "hv")
    title = ggdraw() + draw_label(mains[j], fontface='bold', size = 22)
    g = plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
    print(g)
    
  }
  dev.off()

  return("Completed! Check for four_panels.pdf in your working directory unless otherwise specified in PATH.")
}
