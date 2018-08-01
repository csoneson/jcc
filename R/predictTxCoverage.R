#' Predict coverage profiles of transcripts
#'
#' @param biasModel
#' @param exonsByTx
#' @param genes
#' @param nCores
#' @param verbose
#'
#' @export
#'
#' @author Charlotte Soneson, Michael I Love
#'
#' @importFrom alpine predictCoverage
#' @importFrom GenomicRanges setdiff
#'
predictTxCoverage <- function(biasModel, exonsByTx, bam, genes, nCores,
                              tx2Gene, genome, verbose = FALSE) {

  fitpar <- biasModel
  ebt0 <- exonsByTx

  ## Estimate average fragment length
  avefraglength <- sum(fitpar$`1`$fraglen.density$x * fitpar$`1`$fraglen.density$y/
                         sum(fitpar$`1`$fraglen.density$y))

  ## Load bam file
  bam.files <- bam
  names(bam.files) <- "1"

  ## Get all transcripts
  transcripts <- names(ebt0)
  transcriptsInProvidedGenes <- tx2Gene$tx[tx2Gene$gene %in% genes]
  transcripts <- intersect(transcripts, transcriptsInProvidedGenes)
  names(transcripts) <- transcripts
  if (verbose) message(paste0("Predicting coverage for ", length(transcripts), " transcripts..."))

  res <- mclapply(transcripts, function(tx) {
    message(tx)
    ## Get transcript model
    txmod <- ebt0[[tx]]
    txmods <- sort(ebt0[[tx]])

    ## Predict coverage
    m <- tryCatch({
      a <- alpine::predictCoverage(gene = txmod,
                                   bam.files = bam.files,
                                   fitpar = fitpar,
                                   genome = genome,
                                   model.names = "all")
      if (!is.null(a$`1`$pred.cov$all)) {
        list(covr = a$`1`$pred.cov$all,
             note = "covOK")
      } else {
        ## E.g. if there are no reads mapping to the gene locus
        ## Assume uniform coverage
        list(covr = Rle(rep(1, sum(width(txmod)) - 1)),
             note = "covNA")
      }
    },
    error = function(e) {
      print(paste(c(tx, e), collapse = ": "))
      ## E.g. if the transcript is shorter than the fragment length
      ## Assume uniform coverage
      list(covr = Rle(rep(1, sum(width(txmod)) - 1)),
           note = "covError")
    })

    note <- m$note
    m <- m$covr

    ## Get predicted coverage for junctions
    if (verbose) message("Predicting junction coverages...")
    junctions <- GenomicRanges::setdiff(range(txmods), txmods)
    if (all(strand(txmods) == "+")) {
      junctionpos <- cumsum(width(txmods))
      junctionpos <- junctionpos[-length(junctionpos)]
      mcols(junctions)$pred.cov <- as.numeric(m)[junctionpos]
      strand <- "+"
    } else if (all(strand(txmods) == "-")) {
      junctionpos <- cumsum(width(rev(txmods)))
      junctionpos <- junctionpos[-length(junctionpos)]
      mcols(junctions)$pred.cov <- rev(as.numeric(m)[junctionpos])
      strand <- "-"
    } else {
      strand <- "mixed"
    }

    list(pred.cov = m, strand = strand, junctions = junctions,
         avefraglength = avefraglength, note = note)

  }, mc.preschedule = FALSE, mc.cores = nCores)

  res
}
