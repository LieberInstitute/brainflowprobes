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
#'     \item{Sep }{ phenotype information from
#'     https://www.biorxiv.org/content/10.1101/567966v1; a data.frame with 23
#'     rows and 15 columns. Column descriptions: SampleID is the sample name,
#'     Zone is the RNA fraction of the sample ("Nucleus" or "Cytosol"), Age,
#'     Sex, and Race list these demographic characteristics, Fetal categorizes
#'     each sample age as "Fetal" or "Adult", Library categorizes the RNAseq
#'     library preparation method as polyA selection ("polyA") or rRNA depletion
#'     ("RiboZero"), RIN_fraction is the RNA Integrity Number for each sample,
#'     sumMapped is the number of mapped reads, Shortlabels, Label, col, and
#'     LabelFrac are columns of information for plotting, BigWig is the name of
#'     the BigWig file for each sample, and files lists the URL for the BigWig
#'     online.}
#'     \item{Deg }{ phenotype information from
#'     https://www.pnas.org/content/114/27/7130; a data.frame with 40 rows and
#'     16 columns. Column descriptions: DegradationTime is the number of minutes
#'      the brain tissue for each sample was left on the benchtop at room
#'      temperature before RNA was extracted, AgeDeath is the age of the donor
#'      at time of death, Dx lists whether the donor was diagnosed as
#'      schizophrenic or was a neurotypical control, Sex and Race list these
#'      demographic characteristics for each sample, pH lists the pH of the
#'      brain at collection, PMI is the postmortem interval between death and
#'      brain harvesting, BrNum is the donor number, RIN is the RNA Integrity
#'      Number for the RNA sample, LibraryProtocol categorizes the RNAseq
#'      library preparation method as polyA selection ("polyA") or rRNA
#'      depletion ("Ribo"), SampleID is the ID for the sample, totalAssignedGene
#'       is the proportion of reads mapping to a gene body, sumMapped is the
#'       number of mapped reads, BigWig is the name of the BigWig file,
#'       SampleID_library is the Sample ID and library column values together,
#'       and files lists the URL for the BigWig online.}
#'     \item{Cell }{ phenotype information from
#'     https://www.pnas.org/content/112/23/7285; a data.frame with 466
#'     rows and 11 columns. Column descriptions: geo_accession is the accession
#'     number for each sample in the Gene Expression Omnibus, Age is the numeric
#'      age of each sample, AgeGroup categorizes each sample age as "prenatal"
#'     or "postnatal", sub_tissue lists whether the sample was from cortex or
#'     hippocampus, Cluster_color is the color for each sample for plotting,
#'     Cell_type is the cell identity assigned to each sample, RunName is the
#'     run name assigned to each sample in the Sequence Read Archive (SRA),
#'     SubjectID is the subject label, sumMapped is the number of mapped reads,
#'     BigWig is the name of the BigWig file, and files lists the URL for the
#'     BigWig online.}
#'     \item{Sort }{ phenotype information from
#'     https://www.biorxiv.org/content/10.1101/428391v2; a data.frame with 12
#'     rows and 12 columns. Column descriptions: Description categorizes the
#'     RNAseq library preparation method as polyA selection ("PolyA") or rRNA
#'     depletion ("Ribo"), SubjectID is the subject number, CellType lists
#'     whether the sample was labeled by NeuN antibody (neuronal, "NeuN_Plus")
#'     or not (non-neuronal, "NeuN_Minus"), SampleID combines the sample number,
#'     cell type and library of each sample in one column, BrNum is the donor
#'     ID, RIN is the RNA Integrity Number for each RNA sample, Age is the
#'     numeric age at death, sumMapped is the number of mapped reads, Label and
#'     col provide information for plotting, BigWig is the name of each BigWig
#'     file, and files lists the URL for the BigWig online.}
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

