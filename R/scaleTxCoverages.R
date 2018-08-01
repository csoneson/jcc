#' Scale transcript coverages by abundance estimate
#'
#' @param txCoverageProfiles
#' @param txQuants
#' @param tx2Gene
#' @param strandSpecific
#' @param verbose
#'
#' @export
#'
#' @author Charlotte Soneson
#'
scaleTxCoverages <- function(txCoverageProfiles, txQuants, tx2Gene,
                             strandSpecific, methodName, verbose = FALSE) {

  predcovs <- txCoverageProfiles
  quants <- txQuants
  tx2gene <- tx2Gene

  transcripts <- names(predcovs)
  names(transcripts) <- transcripts

  ## Go through all transcripts and scale predicted junction coverage by the
  ## estimated abundance
  scaledcovs <- lapply(transcripts, function(tx) {
    tryCatch({
      ab <- quants$count[quants$transcript == tx]
      if (length(ab) == 0 || is.na(ab)) ab <- 0  ## if the transcript is not present in the quantification file
      m <- predcovs[[tx]]$junctions
      m$pred.cov <- m$pred.cov / max(1e-10, sum(predcovs[[tx]]$pred.cov)) *
        ab * predcovs[[tx]]$avefraglength
      as.data.frame(m) %>% dplyr::mutate(transcript = tx)
    }, error = function(e) NULL)
  })

  ## Combine estimates for all transcripts
  allcovs <- do.call(rbind, scaledcovs)
  if (!strandSpecific) {
    allcovs$strand <- "*"
  }

  ## Sum coverages of the same junction from different transcripts. Each junction
  ## is present once per gene it is included in.
  allcovs <- allcovs %>%
    dplyr::group_by(seqnames, start, end, strand) %>%
    dplyr::mutate(pred.cov = sum(pred.cov)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(tx2gene %>% dplyr::select(tx, gene, symbol), by = c("transcript" = "tx")) %>%
    dplyr::group_by(seqnames, start, end, strand, gene) %>%
    dplyr::mutate(transcript = paste(transcript, collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(method = methodName)

  ## Get coverage "note" for each transcript
  covnotes <- sapply(transcripts, function(tx) {
    predcovs[[tx]]$note
  })

  ## Add gene info to transcript quant table
  quants <- quants %>% dplyr::left_join(tx2gene %>% dplyr::select(tx, gene, symbol),
                                        by = c("transcript" = "tx")) %>%
    dplyr::mutate(method = methodName) %>%
    dplyr::left_join(data.frame(transcript = names(covnotes),
                                covnote = covnotes,
                                stringsAsFactors = FALSE), by = "transcript") %>%
    dplyr::select(transcript, gene, symbol, count, TPM, method, covnote)

  list(scaledcovs = scaledcovs, allcovs = allcovs, quants = quants)

}
