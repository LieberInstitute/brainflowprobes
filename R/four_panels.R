#' Plot expression coverage in four datasets for candidate probe sequences.
#'
#' \code{four_panels} creates four plots for each candidate probe sequence. The
#' first plot (Separation) shows the adjusted read coverage in cytosolic and
#' nuclear RNA from human postmortem cortex. The second plot (Degradation) shows
#' coverage in human cortical samples exposed to room temperature for 0-60
#' minutes. The third plot (Sorted) shows RNA coverage in nuclei that had been
#' sorted based on reactivity to NeuN-antibody, a neuronal marker. NeuN+ samples
#' are enriched for neurons, and NeuN- samples are enriched for non-neurons. The
#' fourth plot (Single Cells) shows the expression coverage in single cells
#' isolated from human temporal lobe.
#'
#' @param REGION Either a single hg19 genomic sequence including the chromosome,
#'   start, end, and optionally strand separated by colons (e.g.,
#'   "chr20:10199446-10288068:+"), or a string of sequences. Must be character.
#'   Chromosome must be proceeded by "chr".
#' @param PATH This parameter indicates the path where the plot(s) should be
#'   saved. "Default" indicated the file will be saved in the working directory.
#' @param JUNCTIONS A logical value indicating if the candidate probe sequence
#'   spans splice junctions (Default=FALSE).
#' @return \code{four_panels} first annotates the input candidate probe
#'   sequence(s) in REGION using \code{\link[bumphunter]{matchGenes}}, and then
#'   cuts the expression coverage for each sequence from each sample in four
#'   different datasets (see the BrainFlow publication for references) using
#'   \code{\link[derfinder]{getRegionCoverage}}. The coverage is normalized to
#'   the total mapped reads per sample and kilobase width of each probe region
#'   before log2 transformation. The four plots are labeled by the dataset and
#'   the plots are topped by the sequence coordinates, sequence width, and the
#'   name of the nearest gene.
#'
#'   A good candidate probe sequence will have several characteristics. In the
#'   Separation data, the sequence should be relatively highly expressed in
#'   nuclear RNA, at least in your age of interest. The sequence should also
#'   show stable expression over the 60 minutes of room temperature exposure in
#'   the Degradation data. The sequence should also be expressed in the
#'   appropriate NeuN fraction (depending on cell type specificity) in the
#'   Sorted dataset, and also be expressed in the right cell type in the Single
#'   Cell dataset.
#'
#'   \code{four_panels} saves the results as four_panels.pdf in the working
#'   directory unless otherwise specified in PATH.
#'
#'   If JUNCTIONS=TRUE, this means that the candidate probe sequence spans
#'   splice junctions. In this case, the character vector of regions should
#'   represent the coordinates of each exon spanned in the sequence. If
#'   JUNCTIONS=TRUE, \code{four_panels} will sum the coverage of each exon and
#'   plot that value for each dataset instead of creating an independent set of
#'   plots for each exon. This is a way to avoid deflating coverage by including
#'   lowly-expressed intron coverage in the plots.
#' @examples
#' four_panels("chr20:10286777-10288069:+")
#'
#' four_panels(c("chr20:10286777-10288069:+",
#'               "chr18:74690788-74692427:-",
#'               "chr19:49932861-49933829:-"))
#'
#' PENK_exons <- c("chr8:57353587-57354496:-",
#'                 "chr8:57358375-57358515:-",
#'                 "chr8:57358985-57359040:-",
#'                 "chr8:57359128-57359292:-")
#' four_panels(PENK_exons, JUNCTIONS=TRUE, PATH = "/path/to/directory/")
#'
#' \dontrun{
#'
#' four_panels("chr20:10286777-10288069:+", PATH = "/path/to/directory")
#' }


four_panels <- function(REGION, PATH="Default", JUNCTIONS=FALSE) {

  gr <- GenomicRanges::GRanges(REGION)
  nearestAnnotation <- bumphunter::matchGenes(x = gr, subject = genes)

  pd <- list(Sep = pdSep, Deg = pdDeg, Cell = pdCell, Sort = pdSort)

  files <- list(Sep = derfinder::rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Sep',
                                          samplepatt='*.bw$', fileterm=NULL),
                Deg = derfinder::rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Deg',
                                          samplepatt='*.bw$', fileterm=NULL),
                Cell = derfinder::rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Cell',
                                           samplepatt='*.bw$', fileterm=NULL),
                Sort = derfinder::rawFiles(datadir='/dcl01/lieber/ajaffe/Amanda/BigWigs/Sort',
                                           samplepatt='*.bw$', fileterm=NULL))
  files <- mapply(function(bw, mpd) bw[match(mpd$BigWig, names(bw))], files, pd, SIMPLIFY = F)

  regionCov <- mapply(function(bw, mpd) derfinder::getRegionCoverage(
                                           regions = gr, totalMapped = mpd$sumMapped,
                                           files = bw),
                      files, pd, SIMPLIFY = F)

  if (JUNCTIONS!=FALSE) {

    regionCov <- lapply(regionCov, function(x) list(do.call(rbind, x)))
    covMat <- lapply(regionCov, function(x) t(sapply(x, colSums)/100))

    bg <- lapply(covMat, function(x) matrix(sum(width(gr)), nc = ncol(x), nr = nrow(x))/1e3)

    coords <- paste0(seqnames(gr)[1], ":", min(start(gr)),"-", max(end(gr)), ":", strand(gr)[1])
    mains <- paste0(coords, " (", sum(width(gr)), " bp)\n", nearestAnnotation$name)

  } else {

    covMat <- lapply(regionCov, function(x) t(sapply(x, colSums)/100))

    bg <- lapply(covMat, function(x) matrix(width(gr), nc = ncol(x), nr = nrow(x))/1e3)

    coords <- paste0(seqnames(gr), ":", start(gr),"-", end(gr), ":", strand(gr))
    mains <- paste0(coords, " (", width(gr), " bp)\n", nearestAnnotation$name)

  }

  covMat <- mapply(function(cv, b) cv/b, covMat, bg, SIMPLIFY = F)
  covMat <- lapply(covMat, function(x) as.matrix(log2(x + 1)))
  theRanges <- range(do.call(cbind, covMat))

  if (PATH=="Default") {

    pdf("./four_panels.pdf", h = 12, w = 9)

  } else {

    pdf(paste0(PATH, "four_panels.pdf"), h = 12, w = 9)

  }

  for(j in 1:length(regionCov[[1]])) {

    SepDat <- data.frame(Cov = covMat$Sep[j,], LabelFrac = pdSep$LabelFrac,
                         Shortlabels = pdSep$Shortlabels)
    Sep <- ggplot(SepDat, aes(x = LabelFrac, y = Cov)) +
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
            legend.background = element_rect(fill = "transparent",
                                             linetype="solid", colour ="black"),
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      guides(fill=guide_legend(ncol=4)) +
      ggtitle("Separation")


    DegDat <- data.frame(Cov = covMat$Deg[j,], DegradationTime = pdDeg$DegradationTime,
                         LibraryProtocol = pdDeg$LibraryProtocol, BrNum = pdDeg$BrNum)
    Deg <- ggplot(DegDat, aes(x = DegradationTime, y = Cov, fill = BrNum, color = BrNum)) +
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
            legend.background = element_rect(fill = "transparent",
                                             linetype="solid", colour ="black"),
            legend.title = element_blank()) +
      guides(linetype=guide_legend(ncol=2), fill = FALSE, color = FALSE) +
      ggtitle("Degradation")


    SortDat <- data.frame(Cov = covMat$Sort[j,],
                         Label = gsub(":", "\n", levels(factor(pdSort$Label))))
    Sort <- ggplot(SortDat, aes(x = Label, y = Cov)) +
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


    CellDat <- data.frame(Cov = covMat$Cell[j,], Cell_type = pdCell$Cell_type)
    Cell <- ggplot(CellDat, aes(x = Cell_type, y = Cov)) +
      theme_bw() +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(fill = Cell_type), pch=21, color="black") +
      scale_fill_manual(values = c("#FFBB78","#FF9896","#C5B0D5",
                                   "#9467BD","#1F77B4","#2CA02C",
                                   "#D62728","#FF7F0E","#AEC7E8","#98DF8A")) +
      ylim(theRanges) +
      ylab("") + xlab("") +
      theme(text = element_text(size = 20),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      guides(fill=FALSE) +
      ggtitle("Single Cells")

    p <- plot_grid(Sep, Deg, Sort, Cell, ncol = 2, align = "hv")
    title <- ggdraw() + draw_label(mains[j], fontface='bold', size = 22)
    g <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
    print(g)

  }
  dev.off()

  return("Completed! Check for four_panels.pdf in your working
         directory unless otherwise specified in PATH.")

}
