gthr <- function(omega, thr) {
  sapply(omega, function(o) {
    if (is.na(o) || o >= (1 - thr)) 1
    else 0
  })
}

## Define help function for calculating score
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

## Define corresponding help function for calculating scaled coverages
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

  list(junctions = junctionCounts,
       genes = txQuantsGene)

}
