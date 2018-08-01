#' Fit fragment bias model with alpine
#'
#' @param gtf gtf file with genomic features. Preferably in Ensembl format
#' @param bam bam file with read alignments to the genome
#' @param organism The organism
#' @param genome A BSgenome object
#' @param genomeVersion Genome version
#' @param version version
#' @param minLength,maxLength Minimum and maximum length of single-isoform genes
#'   used to fit fragment bias model
#' @param minCount,maxCount Minimum and maximum read coverage of single-isoform
#'   genes used to fit fragment bias model
#' @param subsample Whether to subsample the set of single-isoform genes
#'   satisfying the \code{minLength}, \code{maxLength}, \code{minCount} and
#'   \code{maxCount} criteria before fitting fragment bias model
#' @param nbrSubsample If \code{subsample} is \code{TRUE}, how many genes to
#'   subsample
#' @param verbose Whether to print progress messages
#'
#' @export
#'
#' @author Charlotte Soneson, Michael I Love
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts exonsBy
#' @importFrom ensembldb ensDbFromGtf EnsDb exonsBy
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom Rsamtools countBam BamFile ScanBamParam
#' @importFrom alpine getFragmentWidths buildFragtypes fitBiasModels
#'
fitAlpineBiasModel <- function(gtf, bam, organism, genome, genomeVersion,
                               version, minLength = 600, maxLength = 7000,
                               minCount = 500, maxCount = 10000, subsample = TRUE,
                               nbrSubsample = 200, verbose = FALSE) {

  ## Define the temporary directory where the ensDb object will be saved
  tmpdir <- tempdir()

  ## Create TxDb
  if (verbose) message("Creating TxDb object...")
  txdb <-
    tryCatch({
      ensembldb::ensDbFromGtf(gtf, path = tmpdir, organism = organism,
                              genomeVersion = genomeVersion, version = version)
      ensembldb::EnsDb(paste0(tmpdir, "/", gsub("gtf", "sqlite", basename(gtf))))
    }, error = function(e) {
      GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf",
                                       organism = gsub("_", " ", organism))
    })

  ## Get list of transcripts
  if (verbose) message("Getting list of transcripts...")
  if (class(txdb) == "EnsDb") {
    txdf <- GenomicFeatures::transcripts(txdb, return.type = "DataFrame")  ## data frame format
    txps <- GenomicFeatures::transcripts(txdb)  ## GRanges format
  } else {
    txps <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"), use.names = TRUE)
    txps$gene_id <- unlist(txps$gene_id)
    txdf <- S4Vectors::DataFrame(tx_id = as.character(txps$tx_name),
                                 gene_id = as.character(unlist(txps$gene_id)),
                                 tx_seq_start = as.integer(start(txps)),
                                 tx_seq_end = as.integer(end(txps)),
                                 tx_name = as.character(txps$tx_name))
  }

  ## Get list of exons by transcript
  if (verbose) message("Generating list of exons by transcript")
  if (class(txdb) == "EnsDb") {
    ebt0 <- ensembldb::exonsBy(txdb, by = "tx")
  } else {
    ebt0 <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  }

  ## Select genes with a single isoform
  if (verbose) message("Selecting genes with a single isoform...")
  tab <- table(txdf$gene_id)
  one.iso.genes <- names(tab)[tab == 1]
  if (verbose) message(paste0(length(one.iso.genes), " single-isoform genes found."))

  ## Get transcript names for genes with a single isoform
  one.iso.txs <- txdf$tx_id[txdf$gene_id %in% one.iso.genes]

  ## Extract transcripts from one-isoform genes for use in fitting bias model
  ebt.fit <- ebt0[one.iso.txs]
  ebt.fit <- GenomeInfoDb::keepStandardChromosomes(ebt.fit, pruning.mode = "coarse")

  ## Filter short genes and long genes
  if (verbose) message("Filtering out too short and too long transcripts...")
  min.bp <- minLength
  max.bp <- maxLength
  gene.lengths <- sum(width(ebt.fit))

  ebt.fit <- ebt.fit[gene.lengths > min.bp & gene.lengths < max.bp]
  if (verbose) message(paste0(length(ebt.fit), " transcripts remaining"))

  ## Read bam file
  if (verbose) message("Reading bam file...")
  bam.files <- bam
  names(bam.files) <- "1"

  ## Subset to medium-to-high expressed genes
  if (verbose) message("Subsetting to medium-to-highly expressed genes...")
  txps.fit <- sort(txps[names(ebt.fit)])
  cts <- Rsamtools::countBam(Rsamtools::BamFile(bam.files),
                             param = Rsamtools::ScanBamParam(which = txps.fit))
  S4Vectors::mcols(txps.fit)$cts <- cts$records
  txps.fit <- txps.fit[names(ebt.fit)]
  ebt.fit <- ebt.fit[txps.fit$cts > minCount & txps.fit$cts < maxCount]
  if (verbose) message(paste0(length(ebt.fit), " genes retained."))

  ## Subsample genes to use for the fitting of the bias model
  if (subsample) {
    message("Subsampling genes for fitting bias model...")
    set.seed(1)
    ebt.fit <- ebt.fit[sample(length(ebt.fit), max(nbrSubsample, length(ebt.fit)))]
    if (verbose) message(paste0(length(ebt.fit), " genes retained."))
  }

  ## Get fragment width and read length
  message("Getting fragment widths and read lengths...")
  m <- sort(which(sapply(ebt.fit, length) > 1))
  m0 <- 0
  w <- NULL
  while(is.null(w)) {
    w <- tryCatch({
      m0 <- m0 + 1
      alpine::getFragmentWidths(bam.files, ebt.fit[[m[m0]]])
    }, error = function(e) {
      NULL
    })
  }
  quantiles <- quantile(w, c(.025, .975))
  readlength <- median(getReadLength(bam.files))

  ## Names of genes to retain
  gene.names <- names(ebt.fit)
  names(gene.names) <- gene.names

  ## Build fragment types for selected genes
  if (verbose) message("Building fragment types...")
  fragtypes <- lapply(gene.names, function(gene.name) {
    alpine::buildFragtypes(exons = ebt.fit[[gene.name]],
                           genome = Hsapiens,
                           readlength = readlength,
                           minsize = quantiles[1],
                           maxsize = quantiles[2],
                           gc.str = FALSE)
  })

  ## Define models to fit
  if (verbose) message("Fitting bias model...")
  models <- list(
    "all" = list(
      formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
      ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) + gene",
      offset = c("fraglen", "vlmm")
    )
  )

  ## Fit bias model
  fitpar <- lapply(bam.files, function(bf) {
    alpine::fitBiasModels(genes = ebt.fit,
                          bam.file = bf,
                          fragtypes = fragtypes,
                          genome = Hsapiens,
                          models = models,
                          readlength = readlength,
                          minsize = quantiles[1],
                          maxsize = quantiles[2])
  })

  list(biasModel = fitpar, exonsByTx = ebt0, transcripts = txps)
}
