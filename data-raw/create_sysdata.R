## Create TxDb object for the latest Gencode (v31) version in hg19
library('rtracklayer')
gtf_file <- paste0('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/',
    'release_31/GRCh37_mapping/gencode.v31lift37.annotation.gtf.gz')
gencode_gtf <- import(gtf_file)

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
#     gtf_file,
#     organism = 'Homo sapiens',
#     chrominfo = Seqinfo(genome="hg19")
# )
metadata <- GenomicFeatures:::.prepareGFFMetadata(
    file = gtf_file,
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
## Group by gene (instead of by transcript) such that the TSS displayed
## corresponds to the start of a gene
## Otherwise it could be like
## http://useast.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000132639;r=20:10296230-10307418;t=ENST00000495883
## where chr20:10286777-10288069:+' is closest to a TSS from a lncRNA transcript
## of SNAP25 instead of the coding gene TSS
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


## gs is now available through:
## GenomicState::GenomicStateHub(version = '31', genome = 'hg19', filetype = 'GenomicState')[[1]]

## genes is now available through:
## GenomicState::GenomicStateHub(version = '31', genome = 'hg19', filetype = 'AnnotatedGenes')[[1]]

## Save for later use inside the package
# usethis::use_data(genes, gs, internal = FALSE, overwrite = TRUE)



## Save phenotype tables
pd <- list(Sep = pdSep,
        Deg = pdDeg,
        Cell = pdCell,
        Sort = pdSort)
usethis::use_data(pd, internal = FALSE)

## Original code
# load("/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/Quake_singleCell/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
## Original file location from 2014
## /users/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda
## Likely made using a txdb object that used biomaRt to access Ensembl
# gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
# genes = bumphunter::annotateTranscripts(
# txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
# usethis::use_data(genes, gs, internal = FALSE)

## Coverage data used in the examples
four_panels_example_cov <- brainflowprobes_cov('chr20:10286777-10288069:+')
usethis::use_data(four_panels_example_cov, internal = FALSE)
