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


## Function 1: Print relevant info about the REGION

region_info <- function(REGION, CSV=TRUE, SEQ=TRUE, PATH="Default") {
  
  load("./data/region_info_objects.rda")
  gr = GRanges(REGION)
  nearestAnnotation = matchGenes(x = gr, subject = genes)
  
  if (SEQ==TRUE) {
    df = cbind(as.data.frame(gr), nearestAnnotation, 
               Sequence = as.character(getSeq(Hsapiens, gr)))
  } else {
    df = cbind(as.data.frame(gr), nearestAnnotation)
  }
  
  if (CSV==TRUE) {
    
    if (PATH=="Default") {
      write.csv(df, file="./region_info.csv", quote=FALSE, row.names=FALSE)
    } else {
      write.csv(df, file=paste0(PATH, "region_info.csv"), quote=FALSE, row.names=FALSE)
    }
  }

  return(df)
  return("Completed! If CSV=TRUE, check for region_info.csv in your working directory unless otherwise specified in PATH.")

}
