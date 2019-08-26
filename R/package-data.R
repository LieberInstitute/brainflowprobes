#' GRanges object with the annotated genes for hg19
#'
#' A GRanges object with the annotated genes for hg19 using
#' `bumphunter::annotateTranscripts` on the
#' Gencode v31 annotation lifted over to hg19 grouped by gene
#' instead of by transcript.
#'
#' @name genes
#' @docType data
#' @format A GRanges object with the following annotation columns.
#' \describe{
#'     \item{CSS }{ the coding region start (NA for non-coding genes),}
#'     \item{CSE }{ the coding region end (NA for non-coding genes),}
#'     \item{Tx }{ the Gencode Gene ID,}
#'     \item{Geneid }{ the Ensembl Gene ID,}
#'     \item{Gene }{ the Gene symbol,}
#'     \item{Refseq }{ the RefSeq gene ID,}
#'     \item{Nexons }{ the number of exons,}
#'     \item{Exons }{ the exon coordinates.}
#' }
#' See [annotateTranscripts][bumphunter::annotateTranscripts] for more
#' information.
#'
#' @keywords datasets
#' @seealso [four_panels]
#' <https://github.com/LieberInstitute/brainflowprobes/blob/master/data-raw/create_sysdata.R>
#' <https://www.gencodegenes.org/human/release_31lift37.html>
#'
NULL


#' Genomic state object for hg19
#'
#' The derfinder genomic state object for hg19 created using
#' Gencode v31 lifted over to hg19. Only the `fullGenome` portion of it
#' is saved in this package. See [makeGenomicState][derfinder::makeGenomicState]
#' for more information on how to make this type of object.
#'
#' @name gs
#' @docType data
#' @format A genomic state object for the `fullGenome` for hg19
#' using the Gencode v31 annotation lifted over to hg19. The columns are:
#' \describe{
#'     \item{theRegion }{ type of region of the genome: exon, intron, intergenic
#'     ,}
#'     \item{tx_id }{ IntegerList of the transcripts for the TxDb object
#'     made using Gencode v31 lifted over to hg19 (see links below),}
#'     \item{tx_name }{ the Gencode v31 transcript IDs as a CharacterList,}
#'     \item{gene }{IntegerList of the genes for the TxDb object
#'     made using Gencode v31 lifted over to hg19 (see links below),}
#'     \item{symbol }{ the gene symbols in a CharacterList.}
#' }
#'
#' @keywords datasets
#' @seealso [four_panels]
#' <https://github.com/LieberInstitute/brainflowprobes/blob/master/data-raw/create_sysdata.R>
#' <https://www.gencodegenes.org/human/release_31lift37.html>

NULL


#' A list of phenotype data.frames
#'
#' A list of four phenotype data frames used throughput the package.
#'
#' @name pd
#' @docType data
#' @format A list of four data.frames:
#' \describe{
#'     \item{Sep }{ phenotype information for samples from
#'     <https://www.biorxiv.org/content/10.1101/567966v1>; a data.frame with
#'     23 rows and 15 columns. Column descriptions: SampleID is the sample name,
#'     Zone is the RNA fraction of the sample ("Nucleus" or "Cytosol"), Age,
#'     Sex, and Race list these demographic characteristics, Fetal categorizes
#'     each sample age as "Fetal" or "Adult", Library categorizes the RNAseq
#'     library preparation method as polyA selection ("polyA") or rRNA depletion
#'     ("RiboZero"), RIN_fraction is the RNA Integrity Number for each sample,
#'     sumMapped is the number of mapped reads, Shortlabels, Label, col, and
#'     LabelFrac are columns of information for plotting, BigWig is the name of
#'     the BigWig file for each sample, and files lists the URL for the BigWig
#'     online.}
#'     \item{Deg }{ phenotype information for samples from
#'     <https://www.pnas.org/content/114/27/7130>; a data.frame with 40 rows
#' and 16 columns. Column descriptions: DegradationTime is the number of minutes
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
#'     \item{Cell }{ phenotype information for samples from
#'     <https://www.pnas.org/content/112/23/7285>; a data.frame with 466
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
#'     \item{Sort }{ phenotype information  for samples from
#'     <https://www.biorxiv.org/content/10.1101/428391v2>; a data.frame with
#'     12 rows and 12 columns. Column descriptions: Description categorizes the
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
#' @seealso [four_panels] [brainflowprobes_cov]
#' <https://github.com/LieberInstitute/brainflowprobes/blob/master/data-raw/create_sysdata.R>
#' @examples
#' ##  pd <- list(Sep = pdSep, Deg = pdDeg, Cell = pdCell, Sort = pdSort)
#'
NULL


#' Example base-pair region coverage data
#'
#' A list of base-pair region coverage data.frame lists used for exemplifying
#' the package functionality. This is the data extracted for the example region
#' `chr20:10286777-10288069:+` by looping through the phenotype tables
#' stored in `pd` and using the
#' [getRegionCoverage][derfinder::getRegionCoverage] function. This can be
#' reproduced using the [brainflowprobes_cov] function in this package.
#'
#' @name four_panels_example_cov
#' @docType data
#' @format A list with the base-pair coverage output from
#' [brainflowprobes_cov] with the example region used throughout the
#' package documentation ('chr20:10286777-10288069:+').
#' \describe{
#'     \item{Sep }{ base-pair coverage region data.frame list for pd$Sep,}
#'     \item{Deg }{ base-pair coverage region data.frame list for pd$Deg,}
#'     \item{Cell }{ base-pair coverage region data.frame list for pd$Cell,}
#'     \item{Sort }{ base-pair coverage region data.frame list for pd$Sort.}
#' }
#'
#' @keywords datasets
#' @seealso [four_panels] [plot_coverage]
#' @examples
#' if(FALSE) {
#'     ## Takes about 10 minutes to run!
#'     four_panels_example_cov <- brainflowprobes_cov(
#'         'chr20:10286777-10288069:+'
#'     )
#' }
NULL

