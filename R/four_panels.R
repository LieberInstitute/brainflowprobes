#' Plot expression coverage in four datasets for candidate probe sequences.
#'
#' `four_panels` creates four plots for each candidate probe sequence. The
#' first plot (Separation) shows the adjusted read coverage in cytosolic and
#' nuclear RNA from human postmortem cortex. The second plot (Degradation) shows
#' coverage in human cortical samples exposed to room temperature for 0-60
#' minutes. The third plot (Sorted) shows RNA coverage in nuclei that had been
#' sorted based on reactivity to NeuN-antibody, a neuronal marker. NeuN+ samples
#' are enriched for neurons, and NeuN- samples are enriched for non-neurons. The
#' fourth plot (Single Cells) shows the expression coverage in single cells
#' isolated from human temporal lobe.
#'
#' @param REGION Either a single hg19 genomic sequence including the chromosome,
#'   start, end, and optionally strand separated by colons (e.g.,
#'   'chr20:10199446-10288068:+'), or a string of sequences. Must be character.
#'   Chromosome must be proceeded by 'chr'.
#' @param PDF The name of the PDF file. Defaults to
#' `four_panels.pdf`.
#' @param OUTDIR The default directory where `PDF` will be saved to.
#' @param JUNCTIONS A logical value indicating if the candidate probe sequence
#'   spans splice junctions (Default=FALSE).
#' @param COVERAGE The output of [brainflowprobes_cov] for the input
#' `REGION`. Defaults to `NULL` but it can be pre-computed and saved
#' separately since this is the step that takes the longest to run. Also, this
#' is the only step that depends on rtracklayer's functionality for reading
#' BigWig files which does not run on Windows OS. So it could be run on a
#' non-Windows machine, saved, and then shared with Windows users.
#' @param CODING_ONLY A logical vector of length 1 specifying whether to
#' subset [genes] to only the coding genes. That is, whether to subset
#' [genes] by whether they have a non-NA `CSS` value.
#' @param VERBOSE A logical value indicating whether to print updates from the
#' process of loading the data from the BigWig files.
#' @return `four_panels()` first annotates the input candidate probe
#'   sequence(s) in REGION using [bumphunter::matchGenes()], and then
#'   cuts the expression coverage for each sequence from each sample in four
#'   different datasets (see the BrainFlow publication for references) using
#'   [derfinder::getRegionCoverage()]. The coverage is normalized to
#'   the total mapped reads per sample and kilobase width of each probe region
#'   before log2 transformation. The four plots are labeled by the dataset and
#'   the plots are topped by the sequence coordinates, sequence width, and the
#'   name of the nearest gene.
#'
#'   A good candidate probe sequence will have several characteristics. In the
#'   Separation data, the sequence should be relatively highly expressed in
#'   nuclear RNA, at least in your age of interest. The sequence should also
#'   show stable expression over the 60 minutes of room temperature exposure in
#'   the Degradation data. The sequence should also be expressed in the
#'   appropriate NeuN fraction (depending on cell type specificity) in the
#'   Sorted dataset, and also be expressed in the right cell type in the Single
#'   Cell dataset.
#'
#'   `four_panels()` saves the results as four_panels.pdf in a temporary
#'   directory unless otherwise specified with `OUTDIR`.
#'
#'   `if(JUNCTIONS)`, this means that the candidate probe sequence spans
#'   splice junctions. In this case, the character vector of regions should
#'   represent the coordinates of each exon spanned in the sequence.
#'   `if(JUNCTIONS)`, `four_panels()` will sum the coverage of each exon and
#'   plot that value for each dataset instead of creating an independent set of
#'   plots for each exon. This is a way to avoid deflating coverage by including
#'   lowly-expressed intron coverage in the plots.
#' @examples
#'
#' ## Here we use the pre-saved example coverage data such that this example
#' ## will run fast!
#' four_panels('chr20:10286777-10288069:+',
#'     COVERAGE = four_panels_example_cov)
#'
#' \dontrun{
#' ## Without using COVERAGE, this function reads BigWig files from the web
#' ## using rtracklayer and this functionality is not supported on Windows
#' ## machines.
#' if(.Platform$OS.type != 'windows') {
#'     ## This example takes 10 minutes to run!
#'     four_panels('chr20:10286777-10288069:+')
#' }
#'
#' ## These examples will take several minutes to run depending on your
#' ## internet connection
#' four_panels(c(
#'     'chr20:10286777-10288069:+',
#'     'chr18:74690788-74692427:-',
#'     'chr19:49932861-49933829:-'))
#'
#' PENK_exons <- c(
#'     'chr8:57353587-57354496:-',
#'     'chr8:57358375-57358515:-',
#'     'chr8:57358985-57359040:-',
#'     'chr8:57359128-57359292:-')
#'
#' ## General syntax
#' four_panels(PENK_exons, JUNCTIONS=TRUE,
#'     PDF = 'PDF_file.pdf', OUTDIR = '/path/to/directory/')
#'
#' four_panels('chr20:10286777-10288069:+',
#'     PDF = 'PDF_file.pdf', OUTDIR = '/path/to/directory/')
#'
#'
#' ## Explore the effect of changing CODING_ONLY
#' ## Check how gene name changes in the title of the plot
#' ## (everything else stays the same)
#' cov <- brainflowprobes_cov('chr10:135379301-135379311:+')
#' four_panels('chr10:135379301-135379311:+', COVERAGE = cov)
#' four_panels('chr10:135379301-135379311:+', COVERAGE = cov,
#'     PDF = 'coding_only_four_panels', CODING_ONLY = TRUE)
#' }
#' @export
#' @import GenomicRanges bumphunter ggplot2 derfinder RColorBrewer cowplot
#' @importFrom utils browseURL
#' @importFrom grDevices pdf dev.off
#' @author Amanda J Price


four_panels <- function(REGION,
    PDF = "four_panels.pdf",
    OUTDIR = tempdir(),
    JUNCTIONS = FALSE,
    COVERAGE = NULL,
    CODING_ONLY = FALSE,
    VERBOSE = TRUE) {


    ## For R CMD check
    BrNum <- Cell_type <- Cov <- DegradationTime <- Label <- LabelFrac <-
        LibraryProtocol <- Shortlabels <- NULL

    ## Check the PDF file
    pdf_file <- check_pdf(PDF, OUTDIR)

    ## Define the region(s)
    gr <- GenomicRanges::GRanges(REGION)

    ## Compute the nearest annotation
    nearestAnnotation <- get_nearest_annotation(gr, CODING_ONLY)

    ## Compute or check the coverage
    regionCov <- get_region_cov(REGION, COVERAGE, VERBOSE)

    if(JUNCTIONS) {
        regionCov <- lapply(regionCov, function(x) list(do.call(rbind, x)))
        coords <- paste0(GenomicRanges::seqnames(gr)[1],
            ":", min(GenomicRanges::start(gr)),
            "-", max(GenomicRanges::end(gr)),
            ":", GenomicRanges::strand(gr)[1])

    } else {
        coords <- paste0(GenomicRanges::seqnames(gr),
            ":", GenomicRanges::start(gr),
            "-", GenomicRanges::end(gr),
            ":", GenomicRanges::strand(gr))
    }

    covMat <- lapply(regionCov, function(x) t(sapply(x, colSums)/100))

    bg <- lapply(covMat, function(x) matrix(sum(GenomicRanges::width(gr)),
            ncol = ncol(x),
            nrow = nrow(x))/1000)

    mains <- paste0(coords,
        " (", sum(GenomicRanges::width(gr)),
        " bp)\n",
        nearestAnnotation$name)


    covMat <- mapply(function(cv, b) cv/b, covMat, bg, SIMPLIFY = FALSE)
    covMat <- lapply(covMat, function(x) as.matrix(log2(x + 1)))
    theRanges <- range(do.call(cbind, covMat))

    grDevices::pdf(pdf_file, height = 12, width = 10, useDingbats = FALSE)

    result <- vector(mode = "list", length = length(regionCov[[1]]))
    for (j in seq_len(length(regionCov[[1]]))) {

        SepDat <- data.frame(Cov = covMat$Sep[j,
            ], LabelFrac = brainflowprobes::pd$Sep$LabelFrac,
            Shortlabels = brainflowprobes::pd$Sep$Shortlabels)
        Sep <- ggplot2::ggplot(SepDat,
            aes(x = LabelFrac, y = Cov)) +
            theme_bw() +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(size = 2,
                aes(fill = Shortlabels),
                pch = 21,
                color = "black") +
            scale_fill_manual(values = RColorBrewer::brewer.pal(8,
                "Paired")) +
            labs(fill = "") +
            ylim(theRanges) +
            ylab("Log2(Adj Read/kb)") +
            xlab("") +
            theme(text = element_text(size = 20),
                legend.text = element_text(size = 12),
                legend.position = c(0.5,
                    0.15),
                legend.background = element_rect(fill = "transparent",
                    linetype = "solid",
                    colour = "black"),
                legend.title = element_blank(),
                axis.text.x = element_text(angle = 90,
                    hjust = 1)) +
            guides(fill = guide_legend(ncol = 4)) +
            ggtitle("Separation")


        DegDat <- data.frame(Cov = covMat$Deg[j,
            ], DegradationTime = brainflowprobes::pd$Deg$DegradationTime,
            LibraryProtocol = brainflowprobes::pd$Deg$LibraryProtocol,
            BrNum = brainflowprobes::pd$Deg$BrNum)
        Deg <- ggplot2::ggplot(DegDat,
            aes(x = DegradationTime,
                y = Cov,
                fill = BrNum,
                color = BrNum)) +
            theme_bw() +
            geom_line(aes(linetype = LibraryProtocol),
                lwd = 1.3) +
            geom_point(cex = 2,
                pch = 21,
                aes(fill = BrNum),
                color = "black") +
            scale_color_manual(values = RColorBrewer::brewer.pal(5,
                "Dark2")) +
            scale_fill_manual(values = RColorBrewer::brewer.pal(5,
                "Dark2")) +
            ylim(theRanges) +
            ylab("") +
            xlab("") +
            theme(text = element_text(size = 20),
                legend.text = element_text(size = 12),
                legend.position = c(0.5,
                    0.1),
                legend.background = element_rect(fill = "transparent",
                    linetype = "solid",
                    colour = "black"),
                legend.title = element_blank()) +
            guides(linetype = guide_legend(ncol = 2),
                fill = FALSE,
                color = FALSE) +
            ggtitle("Degradation")


        SortDat <- data.frame(Cov = covMat$Sort[j,
            ], Label = gsub(":",
            "\n", levels(factor(brainflowprobes::pd$Sort$Label))))
        Sort <- ggplot2::ggplot(SortDat,
            aes(x = Label,
                y = Cov)) +
            theme_bw() +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(size = 2,
                aes(fill = Label),
                pch = 21,
                color = "black") +
            scale_fill_manual(values = unique(brainflowprobes::pd$Sort$col)) +
            ylim(theRanges) +
            ylab("Log2(Adj Read/kb)") +
            xlab("") +
            theme(text = element_text(size = 20),
                axis.text.x = element_text(angle = 90,
                    vjust = 0.5,
                    hjust = 1)) +
            guides(fill = FALSE) +
            ggtitle("Sorted")


        CellDat <- data.frame(Cov = covMat$Cell[j,
            ], Cell_type = brainflowprobes::pd$Cell$Cell_type)
        Cell <- ggplot2::ggplot(CellDat,
            aes(x = Cell_type,
                y = Cov)) +
            theme_bw() +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(aes(fill = Cell_type),
                pch = 21,
                color = "black") +
            scale_fill_manual(values = c(
                "#FFBB78",
                "#FF9896",
                "#C5B0D5",
                "#9467BD",
                "#1F77B4",
                "#2CA02C",
                "#D62728",
                "#FF7F0E",
                "#AEC7E8",
                "#98DF8A")) +
            ylim(theRanges) +
            ylab("") +
            xlab("") +
            theme(text = element_text(size = 20),
                axis.text.x = element_text(angle = 90,
                    vjust = 0.5,
                    hjust = 1)) +
            guides(fill = FALSE) +
            ggtitle("Single Cells")

        p <- cowplot::plot_grid(Sep, Deg, Sort, Cell, ncol = 2, align = "hv")
        title <- cowplot::ggdraw() +
            cowplot::draw_label(mains[j],
                fontface = "bold",
                size = 22)
        g <- cowplot::plot_grid(title,
            p, ncol = 1,
            rel_heights = c(0.1,
                1))
        print(g)

        result[[j]] <- g
    }
    grDevices::dev.off()

    .view_pdf(pdf_file)
    return(invisible(j))

}
