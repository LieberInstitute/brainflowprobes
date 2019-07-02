## Make genes object


## Make gs object


load("/dcl01/lieber/ajaffe/Amanda/sorted_nuclear_RNA/Quake_singleCell/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")

gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
genes = bumphunter::annotateTranscripts(
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)

devtools::use_data(genes, gs, pdCell, pdDeg, pdSep, pdSort, internal = TRUE)
