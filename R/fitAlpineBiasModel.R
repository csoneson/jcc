#' Fit fragment bias model with alpine
#'
#' This function provides a wrapper around some of the functions from the
#' \code{alpine} package. Given a gtf file and a bam file with reads aligned to
#' the genome, it will find single-isoform genes (with lengths and expression
#' levels within given ranges) and use the observed read coverages to fit a
#' fragment bias model.
#'
#' @param gtf Path to gtf file with genomic features. Preferably in Ensembl
#'   format.
#' @param bam Path to bam file with read alignments to the genome.
#' @param genome A \code{BSgenome} object.
#' @param organism The organism (e.g., 'Homo_sapiens'). This argument will be
#'   passed to \code{ensembldb::ensDbFromGtf}.
#' @param genomeVersion Genome version (e.g., 'GRCh38'). This argument will be
#'   passed to \code{ensembldb::ensDbFromGtf}.
#' @param version The version of the reference annotation (e.g., 90). This
#'   argument will be passed to \code{ensembldb::ensDbFromGtf}.
#' @param minLength,maxLength Minimum and maximum length of single-isoform genes
#'   used to fit fragment bias model.
#' @param minCount,maxCount Minimum and maximum read coverage of single-isoform
#'   genes used to fit fragment bias model.
#' @param subsample Whether to subsample the set of single-isoform genes
#'   satisfying the \code{minLength}, \code{maxLength}, \code{minCount} and
#'   \code{maxCount} criteria before fitting the fragment bias model.
#' @param nbrSubsample If \code{subsample} is \code{TRUE}, the number of genes
#'   to subsample.
#' @param seed If \code{subsample} is \code{TRUE}, the random seed to use to
#'   ensure reproducibility.
#' @param minSize,maxSize Smallest and largest fragment size to consider. One or
#'   both of these can be \code{NULL}, in which case it is estimated as the 2.5
#'   or 97.5 percentile, respectively, of estimated fragment sizes in the
#'   provided data.
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A list with three elements: \describe{ \item{\code{biasModel}:}{The
#'   fitted fragment bias model.} \item{\code{exonsByTx}:}{A \code{GRangesList}
#'   object with exons grouped by transcript.} \item{\code{transcripts}:}{A
#'   \code{GRanges} object with all the reference transcripts.} }
#'
#' @author Charlotte Soneson, Michael I Love
#'
#' @export
#'
#' @references
#' Soneson C, Love MI, Patro R, Hussain S, Malhotra D, Robinson MD:
#' A junction coverage compatibility score to quantify the reliability of
#' transcript abundance estimates and annotation catalogs. bioRxiv
#' doi:10.1101/378539 (2018).
#'
#' Love MI, Hogenesch JB, Irizarry RA: Modeling of RNA-seq fragment sequence
#' bias reduces systematic errors in transcript abundance estimation. Nature
#' Biotechnology 34(12):1287-1291 (2016).
#'
#' @examples
#' \dontrun{
#' gtf <- system.file("extdata/Homo_sapiens.GRCh38.90.chr22.gtf.gz",
#'                    package = "jcc")
#' bam <- system.file("extdata/reads.chr22.bam", package = "jcc")
#' biasMod <- fitAlpineBiasModel(gtf = gtf, bam = bam,
#'                               organism = "Homo_sapiens",
#'                               genome = Hsapiens, genomeVersion = "GRCh38",
#'                               version = 90, minLength = 230,
#'                               maxLength = 7000, minCount = 10,
#'                               maxCount = 10000, subsample = TRUE,
#'                               nbrSubsample = 30, seed = 1, minSize = NULL,
#'                               maxSize = 220, verbose = TRUE)
#' }
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts exonsBy
#' @importFrom ensembldb ensDbFromGtf EnsDb exonsBy
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom Rsamtools countBam BamFile ScanBamParam
#' @importFrom alpine getFragmentWidths buildFragtypes fitBiasModels
#'   getReadLength
#' @importFrom stats quantile median
#'
fitAlpineBiasModel <- function(gtf, bam, organism, genome, genomeVersion,
                               version, minLength = 600, maxLength = 7000,
                               minCount = 500, maxCount = 10000,
                               subsample = TRUE, nbrSubsample = 200,
                               seed = 1, minSize = NULL,
                               maxSize = NULL, verbose = FALSE) {

  ## Define the temporary directory where the ensDb object will be saved
  tmpdir <- tempdir()

  ## Create TxDb
  if (verbose) message("Creating TxDb object...")
  txdb <-
    tryCatch({
      ensembldb::ensDbFromGtf(gtf, path = tmpdir, organism = organism,
                              genomeVersion = genomeVersion, version = version)
      ensembldb::EnsDb(paste0(tmpdir, "/", gsub("gtf", "sqlite",
                                                basename(gtf))))
    }, error = function(e) {
      GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf",
                                       organism = gsub("_", " ", organism))
    })

  ## Get list of transcripts
  if (verbose) message("Getting list of transcripts...")
  if (class(txdb) == "EnsDb") {
    ## data frame format
    txdf <- GenomicFeatures::transcripts(txdb, return.type = "DataFrame")
    ## GRanges format
    txps <- GenomicFeatures::transcripts(txdb)
  } else {
    txps <- GenomicFeatures::transcripts(txdb,
                                         columns = c("tx_name", "gene_id"),
                                         use.names = TRUE)
    txps$gene_id <- unlist(txps$gene_id)
    txdf <- S4Vectors::DataFrame(tx_id = as.character(txps$tx_name),
                                 gene_id = as.character(unlist(txps$gene_id)),
                                 tx_seq_start = as.integer(start(txps)),
                                 tx_seq_end = as.integer(end(txps)),
                                 tx_name = as.character(txps$tx_name))
  }

  ## Get list of exons by transcript
  if (verbose) message("Generating list of exons by transcript...")
  if (class(txdb) == "EnsDb") {
    ebt0 <- ensembldb::exonsBy(txdb, by = "tx")
  } else {
    ebt0 <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  }

  ## Select genes with a single isoform
  if (verbose) message("Selecting genes with a single isoform...")
  tab <- table(txdf$gene_id)
  oneIsoGenes <- names(tab)[tab == 1]
  if (verbose) message(paste0(length(oneIsoGenes),
                              " single-isoform genes found."))

  ## Get transcript names for genes with a single isoform
  oneIsoTxs <- txdf$tx_id[txdf$gene_id %in% oneIsoGenes]

  ## Extract transcripts from one-isoform genes for use in fitting bias model
  ebtFit <- ebt0[oneIsoTxs]
  ebtFit <- GenomeInfoDb::keepStandardChromosomes(ebtFit,
                                                  pruning.mode = "coarse")

  ## Filter short genes and long genes
  if (verbose) message("Filtering out too short and too long transcripts...")
  minBp <- minLength
  maxBp <- maxLength
  geneLengths <- sum(width(ebtFit))
  ebtFit <- ebtFit[geneLengths > minBp & geneLengths < maxBp]
  if (verbose) message(paste0(length(ebtFit), " transcripts remaining."))

  ## Read bam file
  if (verbose) message("Reading bam file...")
  bamFiles <- bam
  names(bamFiles) <- "1"

  ## Subset to medium-to-high expressed genes
  if (verbose) message("Subsetting to medium-to-highly expressed genes...")
  txpsFit <- sort(txps[names(ebtFit)])
  cts <- Rsamtools::countBam(Rsamtools::BamFile(bamFiles),
                             param = Rsamtools::ScanBamParam(which = txpsFit))
  S4Vectors::mcols(txpsFit)$cts <- cts$records
  txpsFit <- txpsFit[names(ebtFit)]
  ebtFit <- ebtFit[txpsFit$cts > minCount & txpsFit$cts < maxCount]
  if (verbose) message(paste0(length(ebtFit), " genes retained."))

  ## Subsample genes to use for the fitting of the bias model
  if (subsample) {
    if (verbose) message("Subsampling genes for fitting bias model...")
    set.seed(seed)
    ebtFit <- ebtFit[sample(length(ebtFit), min(nbrSubsample, length(ebtFit)))]
    if (verbose) message(paste0(length(ebtFit), " genes retained."))
  }

  ## Get fragment width and read length
  message("Getting fragment widths and read lengths...")
  m <- sort(which(vapply(ebtFit, length, numeric(1)) > 1))
  m0 <- 0
  w <- NULL
  while(is.null(w)) {
    w <- tryCatch({
      m0 <- m0 + 1
      alpine::getFragmentWidths(bamFiles, ebtFit[[m[m0]]])
    }, error = function(e) {
      NULL
    })
  }
  quantiles <- stats::quantile(w, c(.025, .975))
  if (is.null(minSize)) minSize <- round(quantiles[1])
  if (is.null(maxSize)) maxSize <- round(quantiles[2])
  readLength <- stats::median(alpine::getReadLength(bamFiles))

  ## Names of genes to retain
  geneNames <- names(ebtFit)
  names(geneNames) <- geneNames

  ## Build fragment types for selected genes
  if (verbose) message("Building fragment types...")
  fragtypes <- lapply(geneNames, function(geneName) {
    alpine::buildFragtypes(exons = ebtFit[[geneName]],
                           genome = genome,
                           readlength = readLength,
                           minsize = minSize,
                           maxsize = maxSize,
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
  fitpar <- lapply(bamFiles, function(bf) {
    alpine::fitBiasModels(genes = ebtFit,
                          bam.file = bf,
                          fragtypes = fragtypes,
                          genome = genome,
                          models = models,
                          readlength = readLength,
                          minsize = minSize,
                          maxsize = maxSize)
  })

  list(biasModel = fitpar, exonsByTx = ebt0, transcripts = txps)
}
