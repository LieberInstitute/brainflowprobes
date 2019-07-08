#' Plot coverage in nuclear and cytoplasmic RNA for candidate probe sequences.
#'
#' \code{plot_coverage} plots the coverage of nuclear and cytoplasmic RNA from
#' adult and prenatal human prefrontal cortex across a candidate or series of
#' candidate probe sequences for BrainFlow. If the sequence spans splice
#' junctions, the plot will include the introns. A good candidate sequence will
#' be highly and evenly expressed in nuclear RNA.
#'
#' @param REGION Either a single hg19 genomic sequence including the chromosome,
#'   start, end, and optionally strand separated by colons (e.g.,
#'   "chr20:10199446-10288068:+"), or a string of sequences. Must be character.
#'   Chromosome must be proceeded by "chr".
#' @param PATH This parameter indicates the path where the plot(s) should be
#'   saved. "Default" indicated the file will be saved in the working directory.
#' @return \code{plot_coverage} plots all input sequences using
#'   \code{\link[derfinderPlot]{plotRegionCoverage}}. It returns a plot for each
#'   input candidate sequence listed in REGION. Each plot includes coverage of
#'   the sequence(s) in nuclear (N) and cytoplasmic (C) RNA isolated from adult
#'   (A) and fetal (F) prefrontal cortex, sequenced using two different library
#'   preparation methods. PolyA+ libraries (P) were generated using selection
#'   for polyadenylated transcripts, and RiboZero (R) libraries were generated
#'   using a ribosomal RNA depletion protocol.
#'
#'   Each plot also shows the overlapping genes beneath the coverage, and the
#'   genomic location. The title lists the nearest gene, the position of the
#'   sequence relative to the gene's canonical transcriptional start site (TSS),
#'   and further annotation information as described in the "region" column from
#'   \code{\link[bumphunter]{matchGenes}}.
#'
#'   \code{plot_coverage} saves the results as regionCoverage_fractionedData.pdf
#'   in the working directory unless otherwise specified in PATH.
#' @examples
#' plot_coverage("chr20:10286777-10288069:+")
#'
#' plot_coverage(c("chr20:10286777-10288069:+",
#'                 "chr18:74690788-74692427:-",
#'                 "chr19:49932861-49933829:-"))
#'
#' candidates <- c("chr20:10286777-10288069:+",
#'                 "chr18:74690788-74692427:-",
#'                 "chr19:49932861-49933829:-")
#' plot_coverage(candidates, PATH = "/path/to/directory/")
#'
#' \dontrun{
#'
#' plot_coverage("chr20:10286777-10288069:+", PATH = "/path/to/directory")
#' }
#' @export


plot_coverage <- function(REGION, PATH="Default") {

  gr <- GenomicRanges::GRanges(REGION)
  nearestAnnotation <- bumphunter::matchGenes(x = gr, subject = genes)
  annotatedRegions <- derfinder::annotateRegions(regions = gr,
                                                 genomicState = gs, minoverlap = 1)

  regionCov <- derfinder::getRegionCoverage(regions = gr,
                                            totalMapped = pdSep$sumMapped, files = pdSep$files)

  if (PATH=="Default") {

    pdf("./regionCoverage_fractionedData.pdf", h = 8, w = 8)

  } else {

    pdf(paste0(PATH, "regionCoverage_fractionedData.pdf"), h = 8, w = 8)

  }

    p <- derfinderPlot::plotRegionCoverage(regions = gr,
                                           regionCoverage = regionCov,
                                           groupInfo = pdSep$Shortlabels,
                                           colors = brewer.pal(8,"Paired"),
                                           nearestAnnotation = nearestAnnotation,
                                           annotatedRegions = annotatedRegions,
                                           whichRegions = 1:length(gr),
                                           ask = FALSE, verbose = FALSE)
    print(p)
    dev.off()

    return("Completed! Check for regionCoverage_fractionedData.pdf
           in your working directory unless otherwise specified in PATH.")

}
