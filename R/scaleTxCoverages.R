#' Scale transcript coverages by abundance estimate
#'
#' This function takes a list of predicted transcript coverage profiles (from
#' \code{\link{predictTxCoverage}}) and corresponding transcript abundance
#' estimates (from, e.g., Salmon, kallisto, StringTie, RSEM, or any other
#' transcript abundance estimation software), and scales the coverage profiles
#' by the abundance estimates in order to determine the predicted number of
#' reads covering each position in the transcript. In addition, it extracts the
#' predicted number of reads spanning each junction in each transcript, and sums
#' up these for each annotated junction, across all transcripts containing the
#' junction.
#'
#' @param txCoverageProfiles A \code{list} of predicted transcript coverage
#'   profiles (output from \code{\link{predictTxCoverage}}).
#' @param txQuants A \code{data.frame} with estimated transcript abundances.
#'   Must contain columns named \code{transcript}, \code{TPM} and \code{count},
#'   all other columns will be ignored.
#' @param tx2gene A \code{data.frame} with transcript-to-gene mapping. Must have
#'   at least two columns, named \code{tx} and \code{gene}.
#' @param strandSpecific Logical, whether or not the data is strand-specific.
#'   Determines whether or not the strand of a junction should be taken into
#'   account when the coverages are summed up across transcripts.
#' @param methodName Name to use for the abundance estimation method in the
#'   output table.
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A list with two elements:
#' \describe{
#' \item{\code{junctionCov}:}{A \code{tibble} with the predicted coverage for
#' each junction, scaled by the estimated transcript abundances and summarized
#' across all transcripts.}
#' \item{\code{txQuants}:}{A \code{tibble} with estimated transcript abundances,
#' including information about whether or not the coverage profile could be
#' predicted ('covOK') or if a uniform coverage was imposed ('covError' or
#' 'covNA').}
#' }
#'
#' @author Charlotte Soneson
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
#' txQuants <- readRDS(system.file("extdata/quant.sub.rds", package = "jcc"))
#' txsc <- scaleTxCoverages(txCoverageProfiles = preds,
#'                          txQuants = txQuants, tx2Gene = tx2gene,
#'                          strandSpecific = TRUE, methodName = "Salmon",
#'                          verbose = TRUE)
#' }
#'
#' @importFrom dplyr mutate group_by ungroup left_join select distinct select rename
#'
scaleTxCoverages <- function(txCoverageProfiles, txQuants, tx2gene,
                             strandSpecific, methodName, verbose = FALSE) {

  ## Keep only transcript ID, count and TPM columns in quantification table
  txQuants <- txQuants %>% dplyr::select(transcript, count, TPM)

  transcripts <- names(txCoverageProfiles)
  names(transcripts) <- transcripts

  ## Go through all transcripts and scale predicted junction coverage by the
  ## estimated abundance
  scaledcovs <- lapply(transcripts, function(tx) {
    tryCatch({
      ab <- txQuants$count[txQuants$transcript == tx]
      if (length(ab) == 0 || is.na(ab)) ab <- 0  ## if the transcript is not present in the quantification file
      m <- txCoverageProfiles[[tx]]$junctions
      m$pred.cov <- m$pred.cov / max(1e-10, sum(txCoverageProfiles[[tx]]$pred.cov)) *
        ab * txCoverageProfiles[[tx]]$avefraglength
      as.data.frame(m) %>% dplyr::mutate(transcript = tx)
    }, error = function(e) NULL)
  })

  ## Combine estimates for all transcripts
  junctionCov <- do.call(rbind, scaledcovs)
  if (!strandSpecific) {
    junctionCov$strand <- "*"
  }

  ## Sum coverages of the same junction from different transcripts. Each junction
  ## is present once per gene it is included in.
  junctionCov <- junctionCov %>%
    dplyr::group_by(seqnames, start, end, strand) %>%
    dplyr::mutate(pred.cov = sum(pred.cov)) %>%
    dplyr::rename(predCov = pred.cov) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(tx2gene %>% dplyr::select(tx, gene), by = c("transcript" = "tx")) %>%
    dplyr::group_by(seqnames, start, end, strand, gene) %>%
    dplyr::mutate(transcript = paste(transcript, collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(method = methodName)

  ## Get coverage "note" for each transcript
  covNotes <- sapply(transcripts, function(tx) {
    txCoverageProfiles[[tx]]$note
  })

  ## Add gene info to transcript quant table
  txQuants <- txQuants %>% dplyr::left_join(tx2gene %>% dplyr::select(tx, gene),
                                        by = c("transcript" = "tx")) %>%
    dplyr::mutate(method = methodName) %>%
    dplyr::left_join(data.frame(transcript = names(covNotes),
                                covNote = covNotes,
                                stringsAsFactors = FALSE), by = "transcript") %>%
    dplyr::select(transcript, gene, count, TPM, method, covNote)

  list(junctionCov = junctionCov, txQuants = txQuants)
}
