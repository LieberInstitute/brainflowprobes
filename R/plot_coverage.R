#' Plot coverage in nuclear and cytoplasmic RNA for candidate probe sequences.
#'
#' `plot_coverage` plots the coverage of nuclear and cytoplasmic RNA from
#' adult and prenatal human prefrontal cortex across a candidate or series of
#' candidate probe sequences for BrainFlow. If the sequence spans splice
#' junctions, the plot will include the introns. A good candidate sequence will
#' be highly and evenly expressed in nuclear RNA.
#'
#' @param PDF The name of the PDF file. Defaults to
#' `regionCoverage_fractionedData.pdf`.
#' @inheritParams four_panels
#' @return `plot_coverage` plots all input sequences using
#'   [derfinderPlot::plotRegionCoverage()]. It returns a plot for each
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
#'   [bumphunter::matchGenes()].
#'
#'   `plot_coverage` saves the results as regionCoverage_fractionedData.pdf
#'   in a temporary directory unless otherwise specified with `OUTDIR`.
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
#' plot_coverage(candidates, PDF = 'PDF_file.pdf',
#'     OUTDIR = '/path/to/directory/')
#'
#' plot_coverage('chr20:10286777-10288069:+',
#'     PDF = 'PDF_file.pdf', OUTDIR = '/path/to/directory/')
#'
#' ## Explore the effect of changing CODING_ONLY
#' ## Check how gene name and distance to TSS changes in the title of the plot
#' ## (everything else stays the same)
#' cov <- brainflowprobes_cov('chr10:135379301-135379311:+')
#' plot_coverage('chr10:135379301-135379311:+', COVERAGE = cov)
#' plot_coverage('chr10:135379301-135379311:+', COVERAGE = cov,
#'     PDF = 'coding_only_plot_coverage', CODING_ONLY = TRUE)
#' }
#'
#' @export
#' @import GenomicRanges bumphunter derfinder derfinderPlot RColorBrewer
#' @importFrom utils browseURL
#' @importFrom grDevices pdf dev.off
#' @importFrom GenomicState GenomicStateHub
#' @author Amanda J Price


plot_coverage <- function(REGION,
    PDF = "regionCoverage_fractionedData.pdf",
    OUTDIR = tempdir(),
    COVERAGE = NULL,
    CODING_ONLY = FALSE,
    VERBOSE = TRUE) {

    ## Check the PDF file
    pdf_file <- check_pdf(PDF, OUTDIR)

    ## Define the region(s)
    gr <- GenomicRanges::GRanges(REGION)

    ## Compute the nearest annotation
    nearestAnnotation <- get_nearest_annotation(gr, CODING_ONLY)

    ## Obtain the GenomicState object from AnnotationHub
    gs <- GenomicState::GenomicStateHub(version = '31', genome = 'hg19',
        filetype = 'GenomicState')[[1]]

    ## Annotate the regions
    annotatedRegions <- derfinder::annotateRegions(regions = gr,
        genomicState = gs[['fullGenome']],
        minoverlap = 1)

    ## Compute or check the coverage (only for Sep)
    regionCov <- get_region_cov(REGION, COVERAGE, VERBOSE,
        PD = brainflowprobes::pd['Sep'])$Sep

    grDevices::pdf(pdf_file, height = 8, width = 8, useDingbats = FALSE)

    derfinderPlot::plotRegionCoverage(regions = gr,
        regionCoverage = regionCov,
        groupInfo = brainflowprobes::pd$Sep$Shortlabels,
        colors = RColorBrewer::brewer.pal(8, "Paired"),
        nearestAnnotation = nearestAnnotation,
        annotatedRegions = annotatedRegions,
        whichRegions = seq_len(length(gr)),
        ask = FALSE,
        verbose = FALSE
    )
    grDevices::dev.off()

    .view_pdf(pdf_file)
    return(invisible(pdf_file))
}

.view_pdf <- function(pdf_file) {
    message(paste0(Sys.time(), " Completed! Check for ", pdf_file, "."))
    if (interactive())
        utils::browseURL(pdf_file)
}
