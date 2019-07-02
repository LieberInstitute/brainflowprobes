#' Print relevant info about candidate probe sequence.
#'
#' \code{region_info} returns annotation of a single potential probe sequence or
#' list of sequences and, if specified, prints the resuts in a .csv file.
#'
#' @param REGION Either a single hg19 genomic sequence including the chromosome,
#'   start, end, and optionally strand separated by colons (e.g.,
#'   "chr20:10199446-10288068:+"), or a string of sequences to be annotated.
#'   Must be character. Chromosome must be proceeded by "chr".
#' @param CSV A logical value indicating if the results should be exported in a
#'   .csv file.
#' @param SEQ A logical value indicating if the base sequence should be
#'   returned.
#' @param PATH If a .csv file is to be exported, this parameter indicates the
#'   path where the file should be saved. "Default" indicated the file will be
#'   saved in the working directory.
#' @return This function annotates all input sequences using
#'   \code{\link[bumphunter]{matchGenes}}. It returns a data frame where each
#'   row is a genomic sequence specified in REGION. The columns
#'   c("seqnames",	"start", "end",	"width", "strand") list the chromosome,
#'   range, sequence length, and strand of the REGION. The columns c("name",
#'   "annotation", "description", "region", "distance", "subregion",
#'   "insideDistance", "exonnumber", "nexons", "UTR", "geneL", "codingL",
#'   "Geneid", "subjectHits") are described in
#'   \code{\link[bumphunter]{matchGenes}} documentation.
#'
#'   If SEQ=TRUE, a column "Sequence" will be included. This is recommended for
#'   sending the probe sequence to be synthesized.
#'
#'   If CSV=TRUE, a .csv file called region_info.csv will be saved to the
#'   working directory unless otherwise specified in PATH.
#' @examples
#' x <- region_info("chr20:10286777-10288069:+", CSV = FALSE)
#' head(x)
#'
#' y <- region_info(c("chr20:10286777-10288069:+",
#'                    "chr18:74690788-74692427:-",
#'                    "chr19:49932861-49933829:-"),
#'                  CSV = FALSE, SEQ = FALSE)
#' head(y)
#'
#' candidates <- c("chr20:10286777-10288069:+",
#'                 "chr18:74690788-74692427:-",
#'                 "chr19:49932861-49933829:-")
#' region_info(candidates, PATH = "/path/to/directory/")
#'
#' \dontrun{
#' region_info("chr20:10286777-10288069:+", PATH = "/path/to/directory") }


region_info <- function(REGION, CSV=TRUE, SEQ=TRUE, PATH="Default") {

  load("./data/region_info_objects.rda")
  gr = GenomicRanges::GRanges(REGION)
  nearestAnnotation = bumphunter::matchGenes(x = gr, subject = genes)
  nearestAnnotation = nearestAnnotation[,which(colnames(nearestAnnotation)!="strand")]

  if (SEQ==TRUE) {
    df <- cbind(as.data.frame(gr), nearestAnnotation,
               Sequence = as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, gr)))
  } else {
    df <- cbind(as.data.frame(gr), nearestAnnotation)
  }

  if (CSV==TRUE) {

    if (PATH=="Default") {
      write.csv(df, file="./region_info.csv", quote=FALSE, row.names=FALSE)
    } else {
      write.csv(df, file=paste0(PATH, "region_info.csv"), quote=FALSE, row.names=FALSE)
    }
    if (!exists(paste0(PATH, "region_info.csv"))) {
      print("Error: check that the specified path exists. (Are you missing a backslash?)")
    }
  }

  return(df)
  return("Completed! If CSV=TRUE, check for region_info.csv in your working
         directory unless otherwise specified in PATH.")

}
