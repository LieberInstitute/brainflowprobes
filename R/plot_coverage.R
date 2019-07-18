#' Plot coverage in nuclear and cytoplasmic RNA for candidate probe sequences.
#'
#' \code{plot_coverage} plots the coverage of nuclear and cytoplasmic RNA from
#' adult and prenatal human prefrontal cortex across a candidate or series of
#' candidate probe sequences for BrainFlow. If the sequence spans splice
#' junctions, the plot will include the introns. A good candidate sequence will
#' be highly and evenly expressed in nuclear RNA.
#'
#' @param PDF The path and name of the PDF file. Defaults to
#' \code{regionCoverage_fractionedData.pdf}.
#' @inheritParams four_panels
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
#'   and further annotation information as described in the 'region' column from
#'   \code{\link[bumphunter]{matchGenes}}.
#'
#'   \code{plot_coverage} saves the results as regionCoverage_fractionedData.pdf
#'   in the working directory unless otherwise specified in PATH.
#' @examples
#'
#' ## Here we use the pre-saved example coverage data such that this example
#' ## will run fast!
#' plot_coverage('chr20:10286777-10288069:+',
#'     COVERAGE = four_panels_example_cov)
#'
#' ## Without using COVERAGE, this function reads BigWig files from the web
#' ## using rtracklayer and this functionality is not supported on Windows
#' ## machines.
#' if(.Platform$OS.type != 'windows') {
#'     plot_coverage('chr20:10286777-10288069:+',
#'         PDF = 'regionCoverage_fractionedData_fromScratch.pdf')
#' }
#'
#' \dontrun{
#'
#' ## These examples will take a few minutes to run!
#' plot_coverage(c('chr20:10286777-10288069:+',
#'                 'chr18:74690788-74692427:-',
#'                 'chr19:49932861-49933829:-'))
#'
#' candidates <- c('chr20:10286777-10288069:+',
#'                 'chr18:74690788-74692427:-',
#'                 'chr19:49932861-49933829:-')
#'
#' ## General syntax:
#' plot_coverage(candidates, PDF = '/path/to/directory/PDF_file.pdf')
#'
#' plot_coverage('chr20:10286777-10288069:+',
#'     PDF = '/path/to/directory/PDF_file.pdf')
#'
#' ## Explore the effect of changing CODING_ONLY
#' ## Check how gene name and distance to TSS changes in the title of the plot
#' ## (everything else stays the same)
#' cov <- brainflowprobes_cov('chr10:135379301-135379311:+')
#' plot_coverage('chr10:135379301-135379311:+', COVERAGE = cov)
#' plot_coverage('chr10:135379301-135379311:+', COVERAGE = cov,
#'     PDF = 'coding_only_plot_coverage', CODING_ONLY = TRUE)
#' }
#' @export
#' @import GenomicRanges bumphunter derfinder derfinderPlot RColorBrewer
#' @importFrom utils browseURL
#' @importFrom grDevices pdf dev.off
#' @author Amanda J Price


plot_coverage <- function(REGION,
    PDF = "regionCoverage_fractionedData.pdf",
    COVERAGE = NULL,
    CODING_ONLY = FALSE,
    VERBOSE = TRUE) {

    pdf_file <- PDF
    if (!grepl("pdf$",
        tolower(PDF)))
        pdf_file <- paste0(PDF,
            ".pdf")
    if (file.exists(pdf_file))
        stop(paste("The file",
            pdf_file,
            "already exists! Rename or erase it before proceeding."))

    gr <- GenomicRanges::GRanges(REGION)
    gr_subject <- if(CODING_ONLY) {
        brainflowprobes::genes[!is.na(brainflowprobes::genes$CSS)]
    } else {
        brainflowprobes::genes
    }
    nearestAnnotation <- bumphunter::matchGenes(x = gr,
        subject = gr_subject)
    annotatedRegions <- derfinder::annotateRegions(regions = gr,
        genomicState = brainflowprobes::gs,
        minoverlap = 1)

    if(is.null(COVERAGE)) {
        regionCov <- brainflowprobes_cov(
            REGION = REGION,
            PD = brainflowprobes::pd['Sep'],
            VERBOSE = VERBOSE
        )
    } else {
        stopifnot(is.list(COVERAGE))
        stopifnot('Sep' %in% names(COVERAGE))
        regionCov <- COVERAGE
    }
    regionCov <- regionCov$Sep

    grDevices::pdf(pdf_file, height = 8,
        width = 8, useDingbats = FALSE)

    derfinderPlot::plotRegionCoverage(regions = gr,
        regionCoverage = regionCov,
        groupInfo = brainflowprobes::pd$Sep$Shortlabels,
        colors = RColorBrewer::brewer.pal(8,
            "Paired"),
        nearestAnnotation = nearestAnnotation,
        annotatedRegions = annotatedRegions,
        whichRegions = seq_len(length(gr)),
        ask = FALSE,
        verbose = FALSE)
    grDevices::dev.off()



    message("Completed! Check for ",
        pdf_file, " in your working directory unless otherwise specified.")
    if (interactive())
        utils::browseURL(pdf_file)
    return(invisible(pdf_file))
}
