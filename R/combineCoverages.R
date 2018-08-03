#' Combine observed and predicted coverages
#'
#' This function combines predicted junction coverages (from
#' \code{\link{scaleTxCoverages}}) and observed junction coverages (obtained
#' from aligning the reads to the genome and counting the number of reads
#' spanning each junction). It also summarized transcript abundances on the gene
#' level and adds information about the number/fraction of uniquely mapping
#' reads and multimapping reads spanning each junction.
#'
#' @param junctionCounts A \code{data.frame} with observed junction coverages.
#'   Should have columns \code{seqnames}, \code{start}, \code{end},
#'   \code{strand}, \code{uniqreads} and \code{mmreads}.
#' @param junctionPredCovs A \code{data.frame} with predicted junction
#'   coverages.
#' @param txQuants A \code{data.frame} with estimated transcript abundances.
#'
#' @return A list with three elements: \describe{ \item{\code{junctionCovs}:}{A
#'   \code{GRanges} object with predicted and observed junction coverages.}
#'   \item{\code{txQuants}:}{A \code{tibble} with estimated transcript
#'   abundances.} \item{\code{geneQuants}:}{A \code{tibble} with summarized gene
#'   abundances.} }
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @references
#' Soneson C, Love MI, Patro R, Hussain S, Malhotra D, Robinson MD:
#' A junction coverage compatibility score to quantify the reliability of
#' transcript abundance estimates and annotation catalogs. bioRxiv
#' doi:10.1101/378539 (2018)
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
#' tx2gene <- readRDS(system.file("extdata/tx2gene.sub.rds", package = "jcc"))
#' predCovProfiles <- predictTxCoverage(biasModel = biasMod$biasModel,
#'                                      exonsByTx = biasMod$exonsByTx,
#'                                      bam = bam, tx2gene = tx2gene,
#'                                      genome = Hsapiens,
#'                                      genes = c("ENSG00000070371",
#'                                                "ENSG00000093010"),
#'                                      nCores = 1, verbose = TRUE)
#' txQuants <- readRDS(system.file("extdata/quant.sub.rds", package = "jcc"))
#' txsc <- scaleTxCoverages(txCoverageProfiles = predCovProfiles,
#'                          txQuants = txQuants, tx2gene = tx2gene,
#'                          strandSpecific = TRUE, methodName = "Salmon",
#'                          verbose = TRUE)
#' jcov <- read.delim(system.file("extdata/sub.SJ.out.tab", package = "jcc"),
#'                    header = FALSE, as.is = TRUE) %>%
#'  setNames(c("seqnames", "start", "end", "strand", "motif", "annot",
#'            "uniqreads", "mmreads", "maxoverhang")) %>%
#'  dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>%
#'  dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>%
#'  dplyr::select(seqnames, start, end, strand, uniqreads, mmreads) %>%
#'  dplyr::mutate(seqnames = as.character(seqnames))
#' combCov <- combineCoverages(junctionCounts = jcov,
#'                             junctionPredCovs = txsc$junctionPredCovs,
#'                             txQuants = txsc$txQuants)
#' }
#'
#' @importFrom dplyr summarize mutate select ungroup group_by left_join distinct
#'   arrange filter everything
#'
combineCoverages <- function(junctionCounts, junctionPredCovs,
                             txQuants) {

  ## Process observed junction counts
  jcovnostrand <- junctionCounts %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::summarize(uniqreads = sum(uniqreads),
                     mmreads = sum(mmreads)) %>%
    dplyr::mutate(strand = "*") %>% dplyr::ungroup() %>%
    dplyr::select(seqnames, start, end, strand, uniqreads, mmreads) %>%
    as.data.frame()
  jcov <- rbind(junctionCounts, jcovnostrand)

  ## Process predicted junction counts
  junctionPredCovs <- junctionPredCovs %>%
    dplyr::group_by(seqnames, start, end, gene) %>%
    dplyr::mutate(transcript = paste(unique(strsplit(
      paste(unique(transcript), collapse = ","),
      ",")[[1]]), collapse = ",")) %>% ## to make sure we have the same
                                       ## transcript combination for a given
                                       ## junction across all methods
    dplyr::ungroup() %>%
    dplyr::mutate(predCov = replace(predCov, is.na(predCov), 0)) %>%
    dplyr::left_join(jcov, by = c("seqnames", "start", "end", "strand")) %>%
    dplyr::mutate(uniqreads = replace(uniqreads, is.na(uniqreads), 0),
                  mmreads = replace(mmreads, is.na(mmreads), 0)) %>%
    dplyr::mutate(fracunique = uniqreads/(uniqreads + mmreads)) %>%
    dplyr::mutate(fracunique = replace(fracunique, is.na(fracunique), 1)) %>%
    dplyr::ungroup()

  ## Assign a junction ID to each junction
  j0 <- junctionPredCovs %>% dplyr::select(seqnames, start, end, gene) %>%
    dplyr::distinct() %>% dplyr::group_by(gene) %>% dplyr::arrange(start) %>%
    dplyr::mutate(junctionid = paste0("J", seq_len(length(start)))) %>%
    dplyr::ungroup()

  junctionPredCovs <- junctionPredCovs %>% dplyr::left_join(j0) %>%
    dplyr::select(junctionid, dplyr::everything(), transcript)

  geneQuants <- txQuants %>% dplyr::group_by(gene, method) %>%
    dplyr::mutate(TPMrel = TPM/sum(TPM)) %>%
    dplyr::mutate(TPMrel = replace(TPMrel, is.na(TPMrel), 0)) %>%
    dplyr::summarize(count = sum(count),
                     TPM = sum(TPM),
                     nbr_expressed_transcripts = sum(TPM > 0),
                     nbr_expressed_transcripts_5p = sum(TPMrel > 0.05),
                     covOKfraction = length(which(covNote == "covOK"))/
                       length(covNote)) %>%
    dplyr::ungroup()

  totjunctionreads <- junctionPredCovs %>%
    dplyr::filter(method == method[1]) %>%
    dplyr::select(gene, uniqreads, mmreads) %>%
    dplyr::group_by(gene) %>% dplyr::summarize(uniqjuncreads = sum(uniqreads),
                                               mmjuncreads = sum(mmreads))
  geneQuants <- dplyr::left_join(geneQuants, totjunctionreads, by = "gene") %>%
    dplyr::mutate(uniqjuncreads = replace(uniqjuncreads,
                                          is.na(uniqjuncreads), 0)) %>%
    dplyr::mutate(mmjuncreads = replace(mmjuncreads,
                                        is.na(mmjuncreads), 0)) %>%
    dplyr::mutate(uniqjuncfraction =
                    uniqjuncreads/(uniqjuncreads + mmjuncreads)) %>%
    dplyr::mutate(uniqjuncfraction = replace(uniqjuncfraction,
                                             is.na(uniqjuncfraction), 1))

  ## Add number of junctions per gene (from the junction summary). If a gene is
  ## not in the junction summary, it means that it doesn't have any annotated
  ## junctions.
  geneQuants <- geneQuants %>%
    dplyr::left_join(junctionPredCovs %>%
                       dplyr::select(gene, method) %>%
                       dplyr::distinct(),
                     by = c("gene", "method"))

  list(junctionCovs = as.data.frame(junctionPredCovs),
       txQuants = as.data.frame(txQuants),
       geneQuants = as.data.frame(geneQuants))
}
