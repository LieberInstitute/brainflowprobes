#' GRanges object with the annotated genes for hg19
#'
#' A GRanges object with the annotated genes for hg19 using
#' \code{bumphunter::annotateTranscripts} on the
#' TxDb.Hsapiens.UCSC.hg19.knownGene annotation. Takes about 40 seconds to
#' re-compute as measured with \code{system.time()}.
#'
#' @name genes
#' @docType data
#' @format A GRanges object with the following annotation columns.
#' \describe{
#'     \item{CSS }{ the coding region start,}
#'     \item{CSE }{ the coding region end,}
#'     \item{Tx }{ the transcript ID used in TxDb.Hsapiens.UCSC.hg19.knownGene,}
#'     \item{Geneid }{ the Entrez Gene ID,}
#'     \item{Gene }{ the Gene symbol,}
#'     \item{Refseq }{ the RefSeq gene ID,}
#'     \item{Nexons }{ the number of exons,}
#'     \item{Exons }{ the exon coordinates.}
#' }
#' See \link[bumphunter]{annotateTranscripts} for more information.
#'
#' @keywords datasets
#' @seealso \link{four_panels}
#' @references TxDb.Hsapiens.UCSC.hg19.knownGene bumphunter
#' @examples
#'
#' ## system.time(
#' ##     genes <- bumphunter::annotateTranscripts(
#' ##         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#' ##     )
#' ## )
#'
NULL


#' Genomic state object for hg19
#'
#' The derfinder genomic state object for hg19 created using
#' TxDb.Hsapiens.UCSC.hg19.knownGene. Only the \code{fullGenome} portion of it
#' is saved in this package. It takes about 200 seconds to re-make as timed
#' using \code{system.time()}. See \link[derfinder]{makeGenomicState} for
#' more information on how to make this type of object.
#'
#' @name gs
#' @docType data
#' @format A genomic state object for the \code{fullGenome} for hg19
#' using the XX annotation.
#'
#' @keywords datasets
#' @seealso \link{four_panels}
#' @references TODO
#' @examples
#' ## to Amanda: looks like the genomic state object was made using
#' ## Ensembl, unlike the genes info. Ideally they should be from the same
#' ## annotation, so you might want to update the genes info
#'
#' ## system.time(
#' ##     gs <- derfinder::makeGenomicState(
#' ##         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#' ##         chrs = paste0('chr', c(1:22, 'X', 'Y'))
#' ##    )$fullGenome
#' ## )
NULL


#' A list of phenotype data.frames
#'
#' A list of four phenotype data frames used throughput the package.
#'
#' @name pd
#' @docType data
#' @format A list of four data.frames with TODO columns.
#' \describe{
#'     \item{Sep }{ phenotype information from TODO,}
#'     \item{Deg }{ phenotype information from TODO,}
#'     \item{Cell }{ phenotype information from TODO,}
#'     \item{Sort }{ phenotype information from TODO.}
#' }
#'
#' @keywords datasets
#' @seealso \link{four_panels} \link{brainflowprobes_cov}
#' @references TODO
#' @examples
#' ##  pd <- list(Sep = pdSep, Deg = pdDeg, Cell = pdCell, Sort = pdSort)
#'
NULL


#' Example base-pair region coverage data
#'
#' A list of base-pair region coverage matrices used for exemplifying the
#' package functionality. This is the data extracted for the example region
#' \code{chr20:10286777-10288069:+} by looping through the phenotype tables
#' stored in \code{pd} and using the \link[derfinder]{getRegionCoverage}
#' function. This can be reproduced using the \link{brainflowprobes_cov}
#' function in this package.
#'
#' @name four_panels_example_cov
#' @docType data
#' @format A list with the base-pair coverage output from
#' \link{brainflowprobes_cov} with the example region used throughout the
#' package documentation ('chr20:10286777-10288069:+').
#' \describe{
#'     \item{Sep }{ base-pair coverage region matrix for pd$Sep,}
#'     \item{Deg }{ base-pair coverage region matrix for pd$Deg,}
#'     \item{Cell }{ base-pair coverage region matrix for pd$Cell,}
#'     \item{Sort }{ base-pair coverage region matrix for pd$Sort.}
#' }
#'
#' @keywords datasets
#' @seealso \link{four_panels} \link{plot_coverage}
#' @examples
#' if(FALSE) {
#'     ## Takes about 10 minutes to run!
#'     four_panels_example_cov <- brainflowprobes_cov(
#'         'chr20:10286777-10288069:+'
#'     )
#' }
NULL

