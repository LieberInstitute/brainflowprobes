
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


## Function 3: Create REGION plot in fractionated data

plot_coverage = function(REGION, PATH="Default") {
  
  load("./data/probe_design_four_plot_phenotype_data.rda")
  load("./data/region_info_objects.rda")
  load("./data/plot_coverage_objects.rda")
  
  gr = GRanges(REGION)
  nearestAnnotation = matchGenes(x = gr, subject = genes)
  annotatedRegions = annotateRegions(regions = gr, genomicState = gs, minoverlap = 1)
  
  files = rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Sep', samplepatt='*.bw$', fileterm=NULL)
  files = files[match(pdSep$BigWig, names(files))]
  regionCov = getRegionCoverage(regions=gr, totalMapped = pdSep$sumMapped, files = files)
  
  if (PATH=="Default") {
    
    pdf("./regionCoverage_fractionedData.pdf", h = 8, w = 8)
    p = plotRegionCoverage(regions = gr, 
                           regionCoverage = regionCov,
                           groupInfo = pdSep$Shortlabels, 
                           colors = brewer.pal(8,"Paired"), 
                           nearestAnnotation = nearestAnnotation,
                           annotatedRegions = annotatedRegions, 
                           whichRegions = 1:length(gr),
                           ask = FALSE, verbose = FALSE)
    print(p)
    dev.off()
  } else {
    pdf(paste0(PATH, "regionCoverage_fractionedData.pdf"), h = 8, w = 8)
    p = plotRegionCoverage(regions = gr, 
                           regionCoverage = regionCov,
                           groupInfo = pdSep$Shortlabels, 
                           colors = brewer.pal(8,"Paired"), 
                           nearestAnnotation = nearestAnnotation,
                           annotatedRegions = annotatedRegions, 
                           whichRegions = 1:length(gr),
                           ask = FALSE, verbose = FALSE)
    print(p)
    dev.off()
  }
  
  return("Completed! Check for regionCoverage_fractionedData.pdf in your working directory unless otherwise specified in PATH.")
}


