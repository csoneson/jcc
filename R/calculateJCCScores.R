#' Weight function
#'
#' Function to determine the weight for a junction in the calculation of the JCC
#' scores and the scaled coverages. This function is used by
#' \code{junctionScore} and \code{scaledCoverage}, and does not have to be
#' called by the user.
#'
#' @param omega Fraction of uniquely mapped reads among those spanning the
#'   junction.
#' @param thr Threshold value. If \code{omega} is above this value, the junction
#'   contributes to the JCC score of the gene. If \code{omega} is below this
#'   value, the junction doesn't contribute.
#'
#' @return Either 1 or 0, representing the weight for the junction.
#'
#' @author Charlotte Soneson
#'
gthr <- function(omega, thr) {
  vapply(omega, function(o) {
    if (is.na(o) || o >= (1 - thr)) 1
    else 0
  }, numeric(1))
}

#' Help function for JCC score calculation
#'
#' Function to calculate the JCC score for a given gene
#'
#' @param uniqreads Vector with the number of uniquely mapping junction-spanning
#'   reads for each junction in the gene.
#' @param mmreads Vector with the number of multimapping junction-spanning reads
#'   for each junction in the gene.
#' @param predcovs Vector with the predicted junction coverages for each
#'   junction in the gene.
#' @param g Junction weight function that takes as input the fraction of
#'   uniquely mapping reads for the junction, and returns a weight for the
#'   junction contribution to the score.
#' @param beta Scaling exponent. If beta=1, predicted coverages are scaled to
#'   have the same (weighted) sum as the observed coverages. If beta=0, no
#'   scaling is done.
#' @param ... Additional parameters passed to \code{g}.
#'
#' @return The value of the JCC score for the gene
#'
#' @author Charlotte Soneson
#'
junctionScore <- function(uniqreads, mmreads, predcovs, g, beta = 1, ...) {
  omega <- uniqreads/(uniqreads + mmreads)
  omega[mmreads == 0] <- 1  ## if there are no multi-mapping reads, all
                            ## reads are considered to be unique
  w1 <- (sum(g(omega, ...) * uniqreads)/sum(g(omega, ...) * predcovs)) ^ beta
  ## w1 can be non-numeric if all g(omega)=0 (not enough uniquely mapping reads
  ## for any junction) or if g(omega)*predCov=0 for all junctions, even if
  ## g(omega)!=0 for some of them (if the predicted coverage is 0 for a junction
  ## that has non-zero uniquely mapping reads). In both these cases, we don't
  ## scale the predicted coverage (i.e., we set w1=1).
  w1[is.na(w1)] <- 1
  w1[!is.finite(w1)] <- 1
  signif(sum(abs(w1 * g(omega, ...) * predcovs - g(omega, ...) *
                   uniqreads))/sum(g(omega, ...) * uniqreads), 2)
}

#' Help function for calculation of scaled junction coverages
#'
#' Function to calculate the scaled junction coverages for a given gene
#'
#' @param uniqreads Vector with the number of uniquely mapping junction-spanning
#'   reads for each junction in the gene.
#' @param mmreads Vector with the number of multimapping junction-spanning reads
#'   for each junction in the gene.
#' @param predcovs Vector with the predicted junction coverages for each
#'   junction in the gene.
#' @param g Junction weight function that takes as input the fraction of
#'   uniquely mapping reads for the junction, and returns a weight for the
#'   junction contribution to the score.
#' @param beta Scaling exponent. If beta=1, predicted coverages are scaled to
#'   have the same (weighted) sum as the observed coverages. If beta=0, no
#'   scaling is done.
#' @param ... Additional parameters passed to \code{g}.
#'
#' @return A vector of scaled junction coverages
#'
#' @author Charlotte Soneson
#'
scaledCoverage <- function(uniqreads, mmreads, predcovs, g, beta = 1, ...) {
  omega <- uniqreads/(uniqreads + mmreads)
  omega[mmreads == 0] <- 1  ## if there are no multi-mapping reads,
                            ## all reads are considered to be unique
  w1 <- (sum(g(omega, ...) * uniqreads)/sum(g(omega, ...) * predcovs)) ^ beta
  ## w1 can be non-numeric if all g(omega)=0 (not enough uniquely mapping reads
  ## for any junction) or if g(omega)*predCov=0 for all junctions, even if
  ## g(omega)!=0 for some of them (if the predicted coverage is 0 for a junction
  ## that has non-zero uniquely mapping reads). In both these cases, we don't
  ## scale the predicted coverage (i.e., we set w1=1).
  w1[is.na(w1)] <- 1
  w1[!is.finite(w1)] <- 1
  w1 * predcovs
}

#' Calculate scaled coverages and JCC scores
#'
#' This function calculates the scaled predicted junction coverage for each
#' junction, and the JCC score for each gene.
#'
#' @param junctionCovs A \code{data.frame} with observed and predicted junction
#'   counts
#' @param geneQuants A \code{data.frame} with gene-level abundances.
#' @param mmFracThreshold Scalar value in [0, 1] indicating the maximal fraction
#'   of multimapping reads a junction can have and still contribute to the
#'   score.
#'
#' @return A list with two elements:
#' \describe{
#' \item{\code{junctionCovs}:}{A \code{data.frame} with predicted and
#' observed junction coverages.}
#' \item{\code{geneScores}:}{A \code{data.frame} with gene scores.}
#' }
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @references
#' Soneson C, Love MI, Patro R, Hussain S, Malhotra D, Robinson MD: A junction
#' coverage compatibility score to quantify the reliability of transcript
#' abundance estimates and annotation catalogs. bioRxiv doi:10.1101/378539
#' (2018)
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
#' jcc <- calculateJCCScores(junctionCovs = combCov$junctionCovs,
#'                           geneQuants = combCov$geneQuants)
#' }
#'
#' @importFrom dplyr group_by mutate ungroup left_join distinct select %>%
#'
calculateJCCScores <- function(junctionCovs, geneQuants,
                               mmFracThreshold = 0.25) {
  ## Calculate score.
  junctionCovs <- junctionCovs %>%
    dplyr::group_by(gene, method) %>%
    dplyr::mutate(jccscore = junctionScore(uniqreads, mmreads, predCov,
                                           g = gthr, beta = 1,
                                           thr = mmFracThreshold)) %>%
    dplyr::mutate(scaled.cov = scaledCoverage(uniqreads, mmreads, predCov,
                                              g = gthr, beta = 1,
                                              thr = mmFracThreshold)) %>%
    dplyr::mutate(methodscore = paste0(method, " (", jccscore, ")")) %>%
    dplyr::ungroup()

  ## Add score to gene table
  geneQuants <- dplyr::left_join(geneQuants,
                                 junctionCovs %>%
                                   dplyr::select(gene, method, jccscore) %>%
                                   dplyr::distinct(),
                                 by = c("gene", "method"))

  list(junctionCovs = junctionCovs, geneScores = geneQuants)
}
