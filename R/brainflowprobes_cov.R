#' Extract coverage data for a set of regions
#'
#' This function extracts the data from the BigWig coverage files that is then
#' used by \link{four_panels}. This function can take a while to run depending
#' on your internet connection. Furthermore, this function relies on
#' functionality in the rtracklayer package for reading BigWig files which
#' does not work in Windows machines. The data extracted by this function is
#' also used by \link{plot_coverage}.
#'
#' @inheritParams four_panels
#' @param PD A list of data.frames with the \code{sumMapped} and \code{files}
#' columns. Defaults to the data included in this package.
#'
#' @return A list with the region coverage matrices used by
#' \link{four_panels} and \link{plot_coverage}.
#' @export
#' @import derfinder GenomicRanges
#' @author Leonardo Collado-Torres
#'
#' @examples
#'
#' ## This function loads data from BigWig files using the rtracklayer package.
#' ## This functionality is not supported on Windows OS machines!
#' if(.Platform$OS.type != 'windows') {
#'
#' ## How long this takes to run will depend on your internet connection.
#' example_cov <- brainflowprobes_cov('chr20:10286777-10288069:+',
#'      PD = lapply(brainflowprobes::pd, head, n = 2)
#' )
#' }
#'
#' ## This is how the example data included in the package was made:
#' \dontrun{
#' ## This can take about 10 minutes to run!
#' four_panels_example_cov <- brainflowprobes_cov('chr20:10286777-10288069:+')
#' }
#'
#' ## If you are interested, you could download all the BigWig files
#' ## in the \code{brainflowprobes::pd} list of data.frames from the
#' ## \code{files} column to your disk. Doing so will greatly increase the
#' ## speed for \code{brainflowprobes_cov} and the functions that depend on
#' ## this data. Then edit \code{brainflowprobes::pd} \code{files} to point to
#' ## your local files.
#'
#' ## Web location of BigWig files
#' lapply(brainflowprobes::pd, function(x) head(x$files))
#'

brainflowprobes_cov <- function(REGION, PD = brainflowprobes::pd,
    VERBOSE = TRUE) {

    stopifnot(all(sapply(PD, is.data.frame)))
    stopifnot(all(sapply(PD, function(x) all(
        c('sumMapped', 'files') %in% colnames(x)))))

    gr <- GenomicRanges::GRanges(REGION)
    regionCov <- lapply(PD,
        function(x) derfinder::getRegionCoverage(regions = gr,
            totalMapped = x$sumMapped,
            files = x$files,
            verbose = VERBOSE))
    return(regionCov)
}
