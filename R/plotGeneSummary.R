#' @importFrom rtracklayer import
#' @importFrom scatterpie geom_scatterpie
#' @importFrom Gviz GeneRegionTrack GenomeAxisTrack DataTrack
#' @importFrom ggplot2 ggplot aes geom_bar theme ylab scale_fill_manual geom_abline xlab theme_bw expand_limits coord_equal facet_wrap guides
#' @importFrom dplyr filter mutate
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges overlapsAny
#' @importFrom grDevices colorRampPalette
#' @importFrom GenomeInfoDb seqnames
#' @importFrom cowplot plot_grid ggdraw draw_image plot_grid get_legend
#'
plotGeneSummary <- function(gtf, junctionCovs, useGene, bwFile,
                            txQuants, txCol = "transcript_id",
                            geneCol = "gene_id", exonCol = "exon_id") {
  # options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)
  options(ucscChromosomeNames = FALSE)

  genemodels <- rtracklayer::import(gtf)
  idx <- match(c(txCol, geneCol, exonCol), colnames(mcols(genemodels)))
  colnames(mcols(genemodels))[idx] <- c("transcript", "gene", "exon")
  S4Vectors::mcols(genemodels)$symbol <- mcols(genemodels)$transcript
  genemodels <- subset(genemodels, type == "exon")

  jl <- junctionCovs %>% dplyr::filter(gene == useGene)

  names(bwFile) <- ""
  bwCond <- structure("g1", names = "")

  ## Gene model track
  if (!("gene_name" %in% colnames(S4Vectors::mcols(genemodels))))
    S4Vectors::mcols(genemodels)$gene_name <- S4Vectors::mcols(genemodels)[[geneCol]]
  gm <- subset(genemodels, tolower(gene) == tolower(useGene) |
                 tolower(gene_name) == tolower(useGene))
  gm <- subset(gm, gene == gene[1])
  id <- unique(gm$gene_name)
  idshow <- paste0(id, " (", unique(gm$gene), ")")
  show_chr <- unique(seqnames(gm))[1]
  gm <- subset(gm, seqnames == show_chr)
  min_coord <- min(start(gm)) - 0.2*(max(end(gm)) - min(start(gm)))
  max_coord <- max(end(gm)) + 0.05*(max(end(gm)) - min(start(gm)))
  gm$transcript <- factor(gm$transcript, levels = unique(gm$transcript))

  ## Other features in the considered region
  gmo <- genemodels[IRanges::overlapsAny(genemodels,
                                         GRanges(seqnames = show_chr,
                                                 ranges = IRanges(start = min_coord,
                                                                  end = max_coord),
                                                 strand = "*"))]
  gmo <- gmo[!(gmo %in% gm)]
  gmo <- GenomicRanges::reduce(gmo)

  ## Define colors
  muted <- c("#DC050C","#E8601C","#7BAFDE","#1965B0","#B17BA6",
             "#882E72","#F1932D","#F6C141","#F7EE55","#4EB265",
             "#90C987","#CAEDAB","#777777")
  txs <- levels(gm$transcript)
  ncols <- nlevels(gm$transcript)
  cols <- grDevices::colorRampPalette(muted)(ncols)
  names(cols) <- txs

  grtr <- Gviz::GeneRegionTrack(gm, showId = TRUE, col = NULL, fill = cols[gm$transcript],
                                name = "", col.title = "black",
                                background.title = "transparent", min.height = 15)
  grtr2 <- Gviz::GeneRegionTrack(gmo, showId = TRUE, col = "black", fill = "white",
                                 name = "", col.title = "black", showId = FALSE,
                                 background.title = "transparent", min.height = 15)

  gtr <- Gviz::GenomeAxisTrack()

  tracks <- c(gtr, grtr, grtr2)

  multiTracks_rnaseq <- lapply(1:length(bwFile), function(i) {
    assign(paste0("rnaseqtr", i),
           Gviz::DataTrack(range = bwFile[i],
                           type = "histogram",
                           chromosome = unique(GenomeInfoDb::seqnames(gm)),
                           col.title = "black",
                           fill = "grey",
                           col = "grey",
                           col.histogram = "grey",
                           fill.histogram = "grey",
                           cex.title = 0,
           ))
  })
  tracks <- c(multiTracks_rnaseq, tracks)

  rn <- round(1e6*runif(1))
  tmpdir <- tempdir()
  png(paste0(tmpdir, "/gviz", rn, ".png"), width = 10.5,
      height = 5.25, unit = "in", res = 400)
  Gviz::plotTracks(tracks, chromosome = show_chr,
                   from = min_coord, to = max_coord, main = idshow,
                   min.width = 0, min.distance = 0, collapse = FALSE)
  dev.off()

  tpms <- ggplot2::ggplot(txQuants %>% dplyr::filter(gene == useGene) %>%
                            dplyr::mutate(
                              transcript = factor(transcript,
                                                  levels = levels(gm$transcript))),
                          ggplot2::aes(x = method, y = TPM, fill = transcript)) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::ylab("Relative TPM") +
    ggplot2::scale_fill_manual(values = cols, name = "") + xlab("") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
                   legend.text = element_text(size = 7))

  for (tt in txs) {
    jl[[tt]] <- grepl(tt, jl$transcript)
  }
  jl[, txs] <- sweep(jl[, txs], 1, rowSums(jl[, txs]), "/")

  jcov <- ggplot2::ggplot() + ggplot2::geom_abline(intercept = 0, slope = 1) +
    scatterpie::geom_scatterpie(aes(x = uniqreads, y = scaledCov, r = max(scaledCov)/13),
                                cols = txs, data = jl, color = NA) +
    ggplot2::facet_wrap(~ methodscore, nrow = 2) +
    ggplot2::coord_equal(ratio = 1) +
    ggplot2::expand_limits(x = range(c(0, jl$scaledCov, jl$uniqreads)),
                           y = range(c(0, jl$scaledCov, jl$uniqreads))) +
    ggplot2::scale_fill_manual(values = cols, name = "") +
    ggplot2::xlab("Number of uniquely mapped reads spanning junction") +
    ggplot2::ylab("Scaled predicted junction coverage") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = element_text(size = 7),
                   legend.text = element_text(size = 7))

  # png(outpng, width = 12, height = 10, unit = "in", res = 400)
  plt <- cowplot::plot_grid(cowplot::ggdraw() +
                              cowplot::draw_image(paste0(tmpdir, "/gviz", rn, ".png")),
                            cowplot::plot_grid(tpms + ggplot2::theme(legend.position = "none"),
                                               jcov + ggplot2::theme(legend.position = "none"),
                                               nrow = 1, labels = c("B", "C"), rel_widths = c(0.9, 1)),
                            cowplot::get_legend(tpms + ggplot2::theme(legend.direction = "horizontal",
                                                                      legend.justification = "center",
                                                                      legend.box.just = "bottom") +
                                                  ggplot2::guides(fill = guide_legend(nrow = length(txs) %/% 8 + 1))),
                            ncol = 1, rel_heights = c(1, 0.8, 0.1),
                            labels = c("A", "", ""))
  # dev.off()

  unlink(paste0(tmpdir, "/gviz", rn, ".png"))
  plt
}
