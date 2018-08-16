context("jcc-calc")

library(dplyr)
library(BSgenome.Hsapiens.NCBI.GRCh38)

gtf <- system.file("extdata/Homo_sapiens.GRCh38.90.chr22.gtf.gz",
                   package = "jcc")
bam <- system.file("extdata/reads.chr22.bam", package = "jcc")
tx2gene <- readRDS(system.file("extdata/tx2gene.sub.rds", package = "jcc"))
txQuants <- readRDS(system.file("extdata/quant.sub.rds", package = "jcc"))

## Fit fragment bias model
biasMod <- fitAlpineBiasModel(gtf = gtf, bam = bam, organism = "Homo_sapiens",
                              genome = Hsapiens, genomeVersion = "GRCh38",
                              version = 90, minLength = 230, maxLength = 7000,
                              minCount = 10, maxCount = 10000,
                              subsample = TRUE, nbrSubsample = 30, seed = 1,
                              minSize = NULL, maxSize = 220,
                              verbose = TRUE)

## Predict transcript coverage profiles
predCovProfiles <- predictTxCoverage(biasModel = biasMod$biasModel,
                                     exonsByTx = biasMod$exonsByTx,
                                     bam = bam, tx2gene = tx2gene,
                                     genome = Hsapiens,
                                     genes = c("ENSG00000070371"),
                                     nCores = 1, verbose = TRUE)

## Scale predicted coverages based on the estimated abundances
txsc <- scaleTxCoverages(txCoverageProfiles = predCovProfiles,
                         txQuants = txQuants, tx2gene = tx2gene,
                         strandSpecific = TRUE, methodName = "Salmon",
                         verbose = TRUE)

jcov <- read.delim(system.file("extdata/sub.SJ.out.tab", package = "jcc"),
                   header = FALSE, as.is = TRUE) %>%
  setNames(c("seqnames", "start", "end", "strand", "motif", "annot",
             "uniqreads", "mmreads", "maxoverhang")) %>%
  dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>%
  dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>%
  dplyr::select(seqnames, start, end, strand, uniqreads, mmreads) %>%
  dplyr::mutate(seqnames = as.character(seqnames))

## Combine the observed and predicted junction coverages
combCov <- combineCoverages(junctionCounts = jcov,
                            junctionPredCovs = txsc$junctionPredCovs,
                            txQuants = txsc$txQuants)

jcc <- calculateJCCScores(junctionCovs = combCov$junctionCovs,
                          geneQuants = combCov$geneQuants)

test_that("Bias model estimation output has the expected form", {
  expect_length(biasMod, 3)
  expect_named(biasMod, c("biasModel", "exonsByTx", "transcripts"))
  expect_is(biasMod$biasModel, "list")
  expect_is(biasMod$exonsByTx, "GRangesList")
  expect_is(biasMod$transcripts, "GRanges")
  expect_equal(names(biasMod$transcripts), biasMod$transcripts$tx_name)
  expect_length(subset(biasMod$transcripts, gene_id == "ENSG00000070371"), 12)
})

test_that("Predicted coverage profile output has the expected form", {
  expect_length(predCovProfiles, 12)
  expect_equal(sort(names(predCovProfiles)),
               sort(subset(biasMod$transcripts, gene_id == "ENSG00000070371")$tx_name))
  expect_length(predCovProfiles[[1]], 5)
  expect_named(predCovProfiles[[1]], c("predCov", "strand", "junctionCov",
                                       "aveFragLength", "note"))
  expect_is(predCovProfiles[[1]], "list")
  expect_is(predCovProfiles[[1]]$predCov, "Rle")
  expect_is(predCovProfiles[[1]]$strand, "character")
  expect_is(predCovProfiles[[1]]$junctionCov, "GRanges")
  expect_is(predCovProfiles[[1]]$aveFragLength, "numeric")
  expect_is(predCovProfiles[[1]]$note, "character")
})

test_that("Scaled junction coverage output has the expected form", {
  expect_length(txsc, 2)
  expect_is(txsc, "list")
  expect_named(txsc, c("junctionPredCovs", "txQuants"))
  expect_length(txsc$txQuants$transcript, 12)
  expect_equal(txsc$txQuants$gene, rep("ENSG00000070371", 12))
  expect_equal(sum(txsc$txQuants$count), 547.001)
  expect_length(txsc$junctionPredCovs$seqnames, 42)
})

test_that("Observed junction coverage has the expected form", {
  expect_is(jcov, "data.frame")
  expect_named(jcov, c("seqnames", "start", "end", "strand",
                       "uniqreads", "mmreads"))
})

test_that("Combined coverage output has the expected form", {
  expect_is(combCov, "list")
  expect_length(combCov, 3)
  expect_named(combCov, c("junctionCovs", "txQuants", "geneQuants"))
  expect_is(combCov$junctionCovs, "data.frame")
  expect_is(combCov$txQuants, "data.frame")
  expect_is(combCov$geneQuants, "data.frame")
  expect_equal(nrow(combCov$txQuants), 12)
  expect_equal(nrow(combCov$geneQuants), 1)
  expect_equal(combCov$geneQuants$count, 547.001)
  expect_equal(sum(combCov$txQuants$count), 547.001)
})

test_that("JCC score calculations are correct", {
  expect_is(jcc, "list")
  expect_length(jcc, 2)
  expect_named(jcc, c("junctionCovs", "geneScores"))
  expect_equal(signif(sum(abs(sum(jcc$junctionCovs$uniqreads)/sum(jcc$junctionCovs$scaledCov) *
                                jcc$junctionCovs$scaledCov - jcc$junctionCovs$uniqreads)) /
                        sum(jcc$junctionCovs$uniqreads), 2), jcc$geneScores$jccscore)
})
