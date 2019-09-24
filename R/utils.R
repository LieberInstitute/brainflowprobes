#' Check the PDF file
#'
#' This function checks that the PDF file does not exist to avoid
#' overwriting a plot the user has previously made.
#'
#' This is an utility function used by \link{four_panels} and
#' \link{plot_coverage}.
#'
#' @inheritParams four_panels
#'
#' @author Leonardo Collado-Torres
#' @return The path to the PDF file if the file doesn't exist. It appends
#' the pdf file extension if it was absent from `PDF`.
#'
#' @export
#' @examples
#'
#' ## Choose a random PDF file
#' PDF <- paste0('test_', stats::runif(1, max = 1e10))
#'
#' ## Initially this works because the output PDF does not exist.
#' PDF <- check_pdf(PDF)
#' ## It also adds the PDF extension if the user didn't supply it.
#' PDF
#'
#' ## Create a dummy PDF file
#' pdf(file = PDF)
#' plot(1, 1)
#' dev.off()
#'
#' ## Now it doesn't work since the PDF file already exists.
#' testthat::expect_error(check_pdf(basename(PDF)),
#'     'already exists! Rename or erase it')
#'

check_pdf <- function(PDF = 'four_panels.pdf', OUTDIR = tempdir()) {
    pdf_file <- file.path(OUTDIR, PDF)
    if (!grepl("pdf$", tolower(pdf_file)))
        pdf_file <- paste0(pdf_file, ".pdf")
    if (file.exists(pdf_file))
        stop(paste("The file",
            pdf_file,
            "\nalready exists! Rename or erase it before proceeding."),
            call. = FALSE)
    return(pdf_file)
}

#' Compute the nearest annotation to the genes in brainflowprobes
#'
#' For a given set of genomic regions, this function computes the nearest
#' annotation information using the \link{genes} object distributed in this
#' package.
#'
#' This is an utility function used by \link{region_info}, \link{four_panels}
#' and \link{plot_coverage}.
#'
#' @param gr A [GenomicRanges::GRanges()][GenomicRanges::GRanges-class] object.
#' @inheritParams four_panels
#'
#' @return The [bumphunter::matchGenes()] output for the annotation information
#' using the `genes` object in this package (subset to only the coding elements
#' if `CODING_ONLY` was set to `TRUE`).
#' @export
#' @author Leonardo Collado-Torres
#' @examples
#'
#' gr <- GenomicRanges::GRanges('chr10:135379301-135379311:+')
#'
#' get_nearest_annotation(gr)
#' get_nearest_annotation(gr, CODING_ONLY = TRUE)
#'
get_nearest_annotation <- function(gr, CODING_ONLY = FALSE) {
    gr_subject <- if(CODING_ONLY) {
        brainflowprobes::genes[!is.na(brainflowprobes::genes$CSS)]
    } else {
        brainflowprobes::genes
    }
    nearestAnnotation <- bumphunter::matchGenes(x = gr, subject = gr_subject)
    return(nearestAnnotation)
}


#' Check or compute the region coverage from the datasets in brainflowprobes
#'
#' This utility function checks the user-provided region data.frame coverage
#' list (`COVERAGE`) or computes a new one using \link{brainflowprobes_cov}.
#' This is used by \link{four_panels} and \link{plot_coverage}.
#'
#' @inheritParams four_panels
#' @inheritParams brainflowprobes_cov
#'
#' @return If `COVERAGE` is provided and all checks pass, then this function
#' returns `COVERAGE`. Otherwise, it computes a new region coverage data.frame
#' list using \link{brainflowprobes_cov}.
#'
#' @export
#' @author Leonardo Collado-Torres
#' @examples
#'
#' ## If all checks pass, then it returns the COVERAGE
#' stopifnot(identical(
#'     get_region_cov(COVERAGE = four_panels_example_cov),
#'     four_panels_example_cov
#' ))
#'
get_region_cov <- function(REGION, COVERAGE = NULL, VERBOSE = TRUE,
    PD = brainflowprobes::pd) {
    if(is.null(COVERAGE)) {
        regionCov <- brainflowprobes_cov(
            REGION = REGION,
            PD = PD,
            VERBOSE = VERBOSE
        )
    } else {
        stopifnot(is.list(COVERAGE))
        if(!all(c('Sep', 'Deg', 'Cell', 'Sort') %in% names(COVERAGE))) {
            stop("'COVERAGE' should be a list with the elements:\n'",
                paste(names(brainflowprobes::pd), collapse = "', '"), "'.",
                call. = FALSE)
        }
        ## Keep the main ones (in case the user shuffled them)
        COVERAGE <- COVERAGE[names(brainflowprobes::pd)]
        if(!all(vapply(COVERAGE, is.list, logical(1)))) {
            stop("Each of the elements of 'COVERAGE' should be a list of\n",
                "region coverage data.frames as created by\n",
                "brainflowprobes_cov().",
                call. = FALSE)
        }
        if(!all(
            vapply(COVERAGE, function(x) is.data.frame(x[[1]]), logical(1)))) {
            stop("Each of the elements of 'COVERAGE' should be a list of\n",
                "region coverage data.frames as created by\n",
                "brainflowprobes_cov().",
                call. = FALSE)
        }
        if(!identical(
            vapply(brainflowprobes::pd, nrow, integer(1)),
            vapply(COVERAGE, function(x) ncol(x[[1]]), integer(1))
        )) {
            stop("Each of the region coverage data.frame lists in 'COVERAGE'\n",
                "should have a column per each of the samples in\n",
                "brainflowprobes::pd.", call. = FALSE)
        }
        if(length(unique(
            vapply(COVERAGE, function(x) nrow(x[[1]]), integer(1)))) != 1) {
            stop("Each of the region coverage data.frames inside 'COVERAGE'\n",
                "should have the same number of rows (1 row per base-pair).",
                call. = FALSE)
        }
        regionCov <- COVERAGE
    }
    return(regionCov)
}
