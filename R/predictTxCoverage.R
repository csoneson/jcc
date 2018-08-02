#' Predict coverage profiles of transcripts
#'
#' This function predicts coverage profiles for each transcript in a specified
#' set of genes (or all transcripts in the reference catalog), based on the
#' fragment bias model estimated by \code{alpine}.
#'
#' @param biasModel Bias model from \code{alpine}. Can be generated using the
#'   \code{\link{fitAlpineBiasModel}} function.
#' @param exonsByTx A \code{GRangesList} grouping exons by transcript. Can be
#'   generated using the \code{\link{fitAlpineBiasModel}} function.
#' @param bam Path to bam file with read alignments to the genome.
#' @param genes Character vector of gene IDs. Coverage profiles will be
#'   estimated for all transcripts annotated to any of these genes. If set to
#'   \code{NULL}, coverage profiles will be estimated for all transcripts in the
#'   reference catalog.
#' @param nCores Integer, number of cores to use for parallel computations.
#' @param tx2gene A \code{data.frame} with transcript-to-gene mapping.
#' @param genome A \code{BSgenome} object.
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A list with one element per transcript. Each element is a list with
#'   five components:
#' \describe{
#' \item{\code{predCov}:}{The predicted coverage profile for the transcript.}
#' \item{\code{strand}:}{The annotated strand for the transcript.}
#' \item{\code{junctionCov}:}{A \code{data.frame} with genomic coordinates and
#' predicted read coverages for each junction in the transcript.}
#' \item{\code{aveFragLength}:}{The estimated average fragment length.}
#' \item{\code{note}:}{Either 'covOK', 'covError' or 'covNA'. If 'covOK',
#' \code{alpine} was able to predict a coverage profile. If 'covError' or
#' 'covNA', the coverage profile could not be predicted, most likely because the
#' transcript is shorter than the fragment length or because no reads in the BAM
#' file overlapped the transcript. In both these cases, we impose a uniform
#' coverage profile for the transcript.}
#' }
#'
#' @author Charlotte Soneson, Michael I Love
#'
#' @export
#'
#' @references
#' Soneson C, Love MI, Patro R, Hussain S, Malhotra D, Robinson MD: A junction coverage compatibility score to quantify the reliability of transcript abundance estimates and annotation catalogs. bioRxiv doi:10.1101/378539 (2018)
#'
#' @examples
#' \dontrun{
#' gtf <- system.file("extdata/Homo_sapiens.GRCh38.90.chr22.gtf.gz",
#' package = "jcc")
#' bam <- system.file("extdata/reads.chr22.bam", package = "jcc")
#' biasMod <- fitAlpineBiasModel(gtf = gtf, bam = bam, organism = "Homo_sapiens",
#'                               genome = Hsapiens, genomeVersion = "GRCh38",
#'                               version = 90, minLength = 230, maxLength = 7000,
#'                               minCount = 10, maxCount = 10000, subsample = TRUE,
#'                               nbrSubsample = 30, seed = 1, minSize = NULL,
#'                               maxSize = 220, verbose = TRUE)
#' tx2gene <- readRDS(system.file("extdata/tx2gene.sub.rds", package = "jcc"))
#' predCovProfiles <- predictTxCoverage(biasModel = biasMod$biasModel,
#'                                      exonsByTx = biasMod$exonsByTx,
#'                                      bam = bam, tx2gene = tx2gene, genome = Hsapiens,
#'                                      genes = c("ENSG00000070371", "ENSG00000244296"),
#'                                      nCores = 1, verbose = TRUE)
#' }
#'
#' @importFrom alpine predictCoverage
#' @importFrom GenomicRanges setdiff
#'
predictTxCoverage <- function(biasModel, exonsByTx, bam, genes, nCores,
                              tx2gene, genome, verbose = FALSE) {

  ## Estimate average fragment length
  aveFragLength <- sum(biasModel$`1`$fraglen.density$x * biasModel$`1`$fraglen.density$y/
                         sum(biasModel$`1`$fraglen.density$y))

  ## Load bam file
  bamFiles <- bam
  names(bamFiles) <- "1"

  ## Get all transcripts
  transcripts <- names(exonsByTx)
  if (!is.null(genes)) {
    transcriptsInProvidedGenes <- tx2gene$tx[tx2gene$gene %in% genes]
    transcripts <- intersect(transcripts, transcriptsInProvidedGenes)
  }
  names(transcripts) <- transcripts
  if (verbose) message(paste0("Predicting coverage for ",
                              length(transcripts), " transcripts..."))

  res <- mclapply(transcripts, function(tx) {
    if (verbose) message(tx)
    ## Get transcript model
    txmod <- exonsByTx[[tx]]
    txmods <- sort(exonsByTx[[tx]])

    ## Predict coverage
    m <- tryCatch({
      a <- alpine::predictCoverage(gene = txmod,
                                   bam.files = bamFiles,
                                   fitpar = biasModel,
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

    list(predCov = m, strand = strand, junctionCov = junctions,
         aveFragLength = aveFragLength, note = note)

  }, mc.preschedule = FALSE, mc.cores = nCores)

  res
}
