## Make genes object


## Make gs object


load("/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/Quake_singleCell/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")

gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
genes = bumphunter::annotateTranscripts(
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)

pd <- list(Sep = pdSep,
        Deg = pdDeg,
        Cell = pdCell,
        Sort = pdSort)

usethis::use_data(genes, gs, pd, internal = FALSE)

four_panels_example_cov <- brainflowprobes_cov('chr20:10286777-10288069:+')
usethis::use_data(four_panels_example_cov, internal = FALSE)
