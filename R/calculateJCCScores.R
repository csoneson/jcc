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
  sapply(omega, function(o) {
    if (is.na(o) || o >= (1 - thr)) 1
    else 0
  })
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
  omega[mmreads == 0] <- 1  ## if there are no multi-mapping reads, all reads are considered to be unique
  w1 <- (sum(g(omega, ...) * uniqreads)/sum(g(omega, ...) * predcovs)) ^ beta
  ## w1 can be non-numeric if all g(omega)=0 (not enough uniquely mapping reads
  ## for any junction) or if g(omega)*pred.cov=0 for all junctions, even if
  ## g(omega)!=0 for some of them (if the predicted coverage is 0 for a junction
  ## that has non-zero uniquely mapping reads). In both these cases, we don't
  ## scale the predicted coverage (i.e., we set w1=1).
  w1[is.na(w1)] <- 1
  w1[!is.finite(w1)] <- 1
  signif(sum(abs(w1 * g(omega, ...) * predcovs - g(omega, ...) * uniqreads))/sum(g(omega, ...) * uniqreads), 2)
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
  omega[mmreads == 0] <- 1  ## if there are no multi-mapping reads, all reads are considered to be unique
  w1 <- (sum(g(omega, ...) * uniqreads)/sum(g(omega, ...) * predcovs)) ^ beta
  ## w1 can be non-numeric if all g(omega)=0 (not enough uniquely mapping reads
  ## for any junction) or if g(omega)*pred.cov=0 for all junctions, even if
  ## g(omega)!=0 for some of them (if the predicted coverage is 0 for a junction
  ## that has non-zero uniquely mapping reads). In both these cases, we don't
  ## scale the predicted coverage (i.e., we set w1=1).
  w1[is.na(w1)] <- 1
  w1[!is.finite(w1)] <- 1
  w1 * predcovs
}

#' Calculate scaled coverages and JCC scores
#'
#' Function to calculate the scaled predicted junction coverage for each
#' junction, and the JCC score for each gene.
#'
#' @param junctionCounts data.frame with observed and predicted junction counts
#' @param txQuantsGene data.frame with gene-level abundances
#' @param mmfracthreshold Scalar value in [0, 1] indicating the maximal fraction
#'   of multimapping reads a junction can have to contribute to the score.
#'
#' @return A list with two elements: (i) \code{junctions}, a data.frame with
#'   observed and predicted junction coverages; (ii) \code{genes}, a data.frame
#'   with gene-level abundances and scores.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @references
#' Soneson C, Love MI, Patro R, Hussain S, Malhotra D, Robinson MD: A junction coverage compatibility score to quantify the reliability of transcript abundance estimates and annotation catalogs. bioRxiv doi:10.1101/378539 (2018)
#'
#' @examples
#'
#' @importFrom dplyr group_by mutate ungroup left_join distinct select
#'
calculateJCCScores <- function(junctionCounts, txQuantsGene, mmfracthreshold = 0.25) {
  ## Calculate score.
  junctionCounts <- junctionCounts %>%
    dplyr::group_by(gene, method) %>%
    dplyr::mutate(score = junctionScore(uniqreads, mmreads, pred.cov, g = gthr,
                                        beta = 1, thr = mmfracthreshold)) %>%
    dplyr::mutate(scaled.cov = scaledCoverage(uniqreads, mmreads, pred.cov, g = gthr,
                                              beta = 1, thr = mmfracthreshold)) %>%
    dplyr::mutate(methodscore = paste0(method, " (", score, ")")) %>%
    dplyr::ungroup()

  ## Add score to gene table
  txQuantsGene <- dplyr::left_join(txQuantsGene,
                                   junctionCounts %>% dplyr::select(gene, method, score) %>%
                                     dplyr::distinct(),
                                   by = c("gene", "method"))

  list(junctions = junctionCounts, genes = txQuantsGene)
}
