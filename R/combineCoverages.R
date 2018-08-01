#' jcov <- read.delim(junctioncovSTAR, header = FALSE, as.is = TRUE)
#' colnames(jcov) <- c("seqnames", "start", "end", "strand", "motif", "annot",
#'                     "uniqreads", "mmreads", "maxoverhang")
#' jcov <- jcov %>% dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>%
#'   dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>%
#'   dplyr::select(seqnames, start, end, strand, uniqreads, mmreads)

combineCoverages <- function(junctionCounts, predictedCoverages,
                             txQuants) {
  jcovnostrand <- junctionCounts %>% group_by(seqnames, start, end) %>%
    dplyr::summarize(uniqreads = sum(uniqreads),
                     mmreads = sum(mmreads)) %>%
    dplyr::mutate(strand = "*") %>% dplyr::ungroup() %>%
    dplyr::select(seqnames, start, end, strand, uniqreads, mmreads) %>%
    as.data.frame()
  jcov <- rbind(junctionCounts, jcovnostrand)

  predictedCoverages <- predictedCoverages %>%
    dplyr::group_by(seqnames, start, end, gene) %>%
    dplyr::mutate(transcript = paste(unique(strsplit(paste(unique(transcript), collapse = ","),
                                                     ",")[[1]]), collapse = ",")) %>% ## to make sure we have the same transcript combination for a given junction across all methods
    dplyr::ungroup() %>%
    dplyr::mutate(pred.cov = replace(pred.cov, is.na(pred.cov), 0)) %>%
    dplyr::left_join(jcov, by = c("seqnames", "start", "end", "strand")) %>%
    dplyr::mutate(uniqreads = replace(uniqreads, is.na(uniqreads), 0),
                  mmreads = replace(mmreads, is.na(mmreads), 0)) %>%
    dplyr::mutate(fracunique = uniqreads/(uniqreads + mmreads)) %>%
    dplyr::mutate(fracunique = replace(fracunique, is.na(fracunique), 1)) %>%
    dplyr::ungroup()

  j0 <- predictedCoverages %>% dplyr::select(seqnames, start, end, gene) %>%
    dplyr::distinct() %>% dplyr::group_by(gene) %>% dplyr::arrange(start) %>%
    dplyr::mutate(junctionid = paste0("J", seq_len(length(start)))) %>%
    dplyr::ungroup()

  predictedCoverages <- predictedCoverages %>% dplyr::left_join(j0) %>%
    dplyr::select(junctionid, everything(), transcript)

  txQuantsGene <- txQuants %>% dplyr::group_by(gene, method) %>%
    dplyr::mutate(TPMrel = TPM/sum(TPM)) %>%
    dplyr::mutate(TPMrel = replace(TPMrel, is.na(TPMrel), 0)) %>%
    dplyr::summarize(count = sum(count),
                     TPM = sum(TPM),
                     nbr_expressed_transcripts = sum(TPM > 0),
                     nbr_expressed_transcripts_5p = sum(TPMrel > 0.05),
                     covOKfraction = length(which(covnote == "covOK"))/length(covnote)) %>%
    dplyr::ungroup()

  totjunctionreads <- predictedCoverages %>% dplyr::filter(method == method[1]) %>%
    dplyr::select(gene, uniqreads, mmreads) %>%
    dplyr::group_by(gene) %>% dplyr::summarize(uniqjuncreads = sum(uniqreads),
                                               mmjuncreads = sum(mmreads))
  txQuantsGene <- dplyr::left_join(txQuantsGene, totjunctionreads, by = "gene") %>%
    dplyr::mutate(uniqjuncreads = replace(uniqjuncreads, is.na(uniqjuncreads), 0)) %>%
    dplyr::mutate(mmjuncreads = replace(mmjuncreads, is.na(mmjuncreads), 0)) %>%
    dplyr::mutate(uniqjuncfraction = uniqjuncreads/(uniqjuncreads + mmjuncreads)) %>%
    dplyr::mutate(uniqjuncfraction = replace(uniqjuncfraction, is.na(uniqjuncfraction), 1))

  ## Add number of junctions per gene (from the junction summary). If a gene is
  ## not in the junction summary, it means that it doesn't have any annotated
  ## junctions.
  txQuantsGene <- txQuantsGene %>%
    dplyr::left_join(predictedCoverages %>%
                       dplyr::select(gene, method) %>%
                       dplyr::distinct(),
                     by = c("gene", "method"))

  list(junctions = predictedCoverages, transcripts = txQuants,
       genes = txQuantsGene)

}
