#' Print relevant info about candidate probe sequence.
#'
#' `region_info` returns annotation of a single potential probe sequence or
#' list of sequences and, if specified, prints the resuts in a .csv file.
#'
#' @param REGION Either a single hg19 genomic sequence including the chromosome,
#'   start, end, and optionally strand separated by colons (e.g.,
#'   `'chr20:10199446-10288068:+'`), or a string of sequences to be annotated.
#'   Must be character. Chromosome must be proceeded by 'chr'.
#' @param CSV A `logical(1)` value indicating if the results should be exported
#' in a .csv file.
#' @param SEQ A `logical(1)`` value indicating if the base sequence should be
#'   returned.
#' @param OUTDIR If a .csv file is to be exported, this parameter indicates the
#'   path where the file should be saved. By default the file will be
#'   saved in a temporary directory.
#' @inheritParams four_panels
#' @return This function annotates all input sequences using
#'   [bumphunter::matchGenes()]. It returns a data frame where each
#'   row is a genomic sequence specified in REGION. The columns
#'   c('seqnames', 'start', 'end', 'width', 'strand') list the chromosome,
#'   range, sequence length, and strand of the REGION. The columns c('name',
#'   'annotation', 'description', 'region', 'distance', 'subregion',
#'   'insideDistance', 'exonnumber', 'nexons', 'UTR', 'geneL', 'codingL',
#'   'Geneid', 'subjectHits') are described in
#'   [bumphunter::matchGenes()] documentation.
#'
#'   If SEQ=TRUE, a column 'Sequence' will be included. This is recommended for
#'   sending the probe sequence to be synthesized.
#'
#'   If CSV=TRUE, a .csv file called region_info.csv will be saved to a
#'   temporary directory unless otherwise specified in `OUTDIR`.
#' @examples
#' x <- region_info('chr20:10286777-10288069:+', CSV = FALSE)
#' head(x)
#'
#' ## You can easily transform this data.frame to a GRanges object
#' GenomicRanges::GRanges(x)
#'
#' y <- region_info(c('chr20:10286777-10288069:+',
#'                    'chr18:74690788-74692427:-',
#'                    'chr19:49932861-49933829:-'),
#'                  CSV = FALSE, SEQ = FALSE)
#' head(y)
#'
#' candidates <- c('chr20:10286777-10288069:+',
#'                 'chr18:74690788-74692427:-',
#'                 'chr19:49932861-49933829:-')
#' region_info(candidates, CSV = FALSE)
#'
#' ## Explore the effect of changing CODING_ONLY
#' ## Check how the "distance", "name", "Geneid" among other values change
#' region_info('chr10:135379301-135379311:+', CSV = FALSE)
#' region_info('chr10:135379301-135379311:+', CSV = FALSE, CODING_ONLY = TRUE)
#'
#'
#' \dontrun{
#' region_info(candidates, OUTDIR = '/path/to/directory/')
#'
#' region_info('chr20:10286777-10288069:+', OUTDIR = '/path/to/directory')
#'
#' }
#' @export
#' @import GenomicRanges bumphunter Biostrings BSgenome.Hsapiens.UCSC.hg19
#' @importFrom utils write.csv
#' @author Amanda J Price


region_info <- function(REGION, CSV = TRUE, SEQ = TRUE, OUTDIR = tempdir(),
    CODING_ONLY = FALSE) {

    ## Define the region(s)
    gr <- GenomicRanges::GRanges(REGION)

    ## Compute the nearest annotation
    nearestAnnotation <- get_nearest_annotation(gr, CODING_ONLY)
    nearestAnnotation <- nearestAnnotation[,
        -which(colnames(nearestAnnotation) %in%
            c("strand", "subjectHits"))]

    if (SEQ) {
        df <- cbind(as.data.frame(gr),
            nearestAnnotation,
            Sequence = as.character(Biostrings::getSeq(
                BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                gr)))
    } else {
        df <- cbind(as.data.frame(gr),
            nearestAnnotation)
    }

    if (CSV) {
        csv_path <- file.path(OUTDIR,
            "region_info.csv")
        if (file.exists(csv_path))
            stop(paste("The file",
                csv_path,
                "already exists! Rename or erase it before proceeding."))
        utils::write.csv(df,
            file = csv_path,
            quote = FALSE,
            row.names = FALSE)

        if (!file.exists(csv_path)) {
            stop(paste("Check that the specified path exists.",
                "(Are you missing a backslash?)"))
        }
    }

    message(
        "Completed! If CSV=TRUE, check for region_info.csv in the temporary\n",
        "directory (i.e. tempdir()) unless otherwise specified in OUTDIR.")
    return(df)

}
