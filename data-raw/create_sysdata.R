## Create TxDb object for the latest Ensembl version in hg19
library('rtracklayer')
ens_gtf <- import('ftp://ftp.ensembl.org/pub/grch37/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz')

library('GenomeInfoDb')
seqlevelsStyle(ens_gtf) <- 'UCSC'

library('GenomicFeatures')
## Doesn't work because of the different seqlevelStyles
# txdb <- makeTxDbFromGFF(
#     'ftp://ftp.ensembl.org/pub/grch37/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz',
#     organism = 'Homo sapiens',
#     chrominfo = Seqinfo(genome="hg19")
# )

metadata <- GenomicFeatures:::.prepareGFFMetadata(
    file = 'ftp://ftp.ensembl.org/pub/grch37/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz',
    dataSource = NA, organism = 'Homo sapiens',
    taxonomyId = NA, miRBaseBuild = NA, metadata = NULL)
gr <- GenomicFeatures:::.tidy_seqinfo(
    gr = ens_gtf,
    circ_seqs = GenomicFeatures::DEFAULT_CIRC_SEQS,
    chrominfo = Seqinfo(genome="hg19")
)
txdb <- makeTxDbFromGRanges(gr, metadata = metadata)



## Extract the genes
library('bumphunter')
genes <- bumphunter::annotateTranscripts(txdb,
    mappingInfo = list('column' = 'ENTREZID', 'keytype' = 'ENSEMBL', 'multiVals' = 'first')
)

## Now build the genomic state object
library('derfinder')
GenomicState.Hsapiens.ensembl.GRCh37.v87 <- makeGenomicState(txdb)
gs <- GenomicState.Hsapiens.ensembl.GRCh37.v87$fullGenome


## Add the symbols
library('org.Hs.eg.db')
gene_gr <- genes(txdb)
gene_gr$symbol <- mapIds(org.Hs.eg.db, gene_gr$gene_id, 'SYMBOL', 'ENSEMBL')
table(is.na(gene_gr$symbol))
# FALSE  TRUE
# 24496 33240
gs$symbol <- extractList(gene_gr$symbol, gs$gene)

## Save for later use inside the package
usethis::use_data(genes, gs, internal = FALSE, overwrite = TRUE)

genes_ens <- genes
gs_ens <- gs




## Create TxDb object for the latest Gencode (v31) version in hg19
library('rtracklayer')
gencode_gtf <- import('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gtf.gz')

library('GenomeInfoDb')
## Keep only the main chrs
gencode_gtf <- keepSeqlevels(
    gencode_gtf,
    paste0('chr', c(1:22, 'X', 'Y', 'M')),
    pruning.mode = 'coarse'
)

library('GenomicFeatures')
# Doesn't work because of the different seqlevels
# txdb <- makeTxDbFromGFF(
#     'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gtf.gz',
#     organism = 'Homo sapiens',
#     chrominfo = Seqinfo(genome="hg19")
# )

metadata <- GenomicFeatures:::.prepareGFFMetadata(
    file = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gtf.gz',
    dataSource = NA, organism = 'Homo sapiens',
    taxonomyId = NA, miRBaseBuild = NA, metadata = NULL)
gr <- GenomicFeatures:::.tidy_seqinfo(
    gr = gencode_gtf,
    circ_seqs = GenomicFeatures::DEFAULT_CIRC_SEQS,
    chrominfo = Seqinfo(genome="hg19")
)
txdb <- makeTxDbFromGRanges(gr, metadata = metadata)



## Extract the genes
library('bumphunter')
genes <- bumphunter::annotateTranscripts(txdb,
    by = 'gene',
    mappingInfo = list('column' = 'ENTREZID', 'keytype' = 'ENSEMBL', 'multiVals' = 'first'),
    simplifyGeneID = TRUE
)

## Now build the genomic state object
library('derfinder')
GenomicState.Hsapiens.gencode.GRCh37.v31 <- makeGenomicState(txdb)
gs <- GenomicState.Hsapiens.gencode.GRCh37.v31$fullGenome


## Add the symbols
library('org.Hs.eg.db')
gene_gr <- genes(txdb)
gene_gr$symbol <- mapIds(
    org.Hs.eg.db,
    gsub('\\..*', '', gene_gr$gene_id),
    'SYMBOL', 'ENSEMBL')
table(is.na(gene_gr$symbol))
# FALSE  TRUE
# 25178 37015
stopifnot(max(unlist(gs$gene)) == length(gene_gr))
gs$symbol <- extractList(gene_gr$symbol, gs$gene)

## Save for later use inside the package
usethis::use_data(genes, gs, internal = FALSE, overwrite = TRUE)




xx <- makeTxDbPackageFromBiomart(
    version = '0.99',
    maintainer = 'Leonardo Collado Torres',
    author='Leonardo Collado Torres',
    biomart = 'ENSEMBL_MART_ENSEMBL',
    dataset = 'hsapiens_gene_ensembl'
)




txdb <- makeTxDbFromEnsembl(organism="Homo sapiens", release = 87)

lookup_dbname_forcehg19 <- function (organism, release = NA)
{
    organism <- GenomicFeatures:::.normarg_organism(organism)
    if (!isSingleNumberOrNA(release))
        stop("'release' must be a valid Ensembl release number or NA")
    available_dbs <- GenomicFeatures:::Ensembl_listMySQLCoreDirs(release = release, use.grch37 = TRUE)
    prefix <- paste0(gsub(" ", "_", tolower(organism), fixed = TRUE),
        "_core_")
    i <- match(prefix, substr(available_dbs, 1L, nchar(prefix)))
    dbname <- available_dbs[[i]]
    if (!is.na(release))
        stopifnot(GenomicFeatures:::.dbname2release(dbname) == release)
    dbname
}


makeTxDbFromEnsembl_forcehg19 <- function (organism = "Homo sapiens", release = NA, circ_seqs = GenomicFeatures::DEFAULT_CIRC_SEQS,
    server = "ensembldb.ensembl.org", username = "anonymous",
    password = NULL, port = 0L, tx_attrib = NULL)
{
    if (!requireNamespace("RMariaDB", quietly = TRUE))
        stop(wmsg("Couldn't load the RMariaDB package. ", "You need to install the RMariaDB package ",
            "in order to use makeTxDbFromEnsembl()."))
    dbname <- lookup_dbname_forcehg19(organism, release = release)
    dbconn <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = dbname,
        host = server, username = username, password = password,
        port = port)
    on.exit(DBI::dbDisconnect(dbconn))
    transcripts <- GenomicFeatures:::.fetch_Ensembl_transcripts(dbconn, tx_attrib)
    splicings <- GenomicFeatures:::.fetch_Ensembl_splicings(dbconn, tx_attrib)
    seq_region_ids <- unique(c(transcripts$seq_region_id, splicings$seq_region_id))
    chrominfo <- GenomicFeatures:::.fetch_Ensembl_chrominfo(dbconn, seq_region_ids = seq_region_ids,
        circ_seqs = circ_seqs)
    m <- match(transcripts$seq_region_id, chrominfo$seq_region_id)
    transcripts$tx_chrom <- chrominfo$chrom[m]
    transcripts$seq_region_id <- NULL
    m <- match(splicings$seq_region_id, chrominfo$seq_region_id)
    splicings$exon_chrom <- chrominfo$chrom[m]
    splicings$seq_region_id <- NULL
    chrominfo$seq_region_id <- NULL
    metadata <- GenomicFeatures:::.gather_Ensembl_metadata(organism, dbname, server)
    message("Make the TxDb object ... ", appendLF = FALSE)
    txdb <- GenomicFeatures::makeTxDb(transcripts, splicings, chrominfo = chrominfo,
        metadata = metadata, reassign.ids = TRUE)
    message("OK")
    txdb
}

txdb <- makeTxDbFromEnsembl_forcehg19(organism="Homo sapiens", release = 97)


library('GenomeInfoDb')
seqlevelsStyle(txdb) <- 'UCSC'

txdb <- loadDb(file.path('TxDb.Hsapiens.BioMart.ensembl.GRCh37.p11', 'inst',
    'extdata', 'TxDb.Hsapiens.BioMart.ensembl.GRCh37.p11.sqlite'))


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
