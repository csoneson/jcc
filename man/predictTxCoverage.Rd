% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictTxCoverage.R
\name{predictTxCoverage}
\alias{predictTxCoverage}
\title{Predict coverage profiles of transcripts}
\usage{
predictTxCoverage(biasModel, exonsByTx, bam, genes, nCores, tx2gene, genome,
  verbose = FALSE)
}
\arguments{
\item{biasModel}{Bias model from \code{alpine}. Can be generated using the
\code{\link{fitAlpineBiasModel}} function.}

\item{exonsByTx}{A \code{GRangesList} grouping exons by transcript. Can be
generated using the \code{\link{fitAlpineBiasModel}} function.}

\item{bam}{Path to bam file with read alignments to the genome.}

\item{genes}{Character vector of gene IDs. Coverage profiles will be
estimated for all transcripts annotated to any of these genes. If set to
\code{NULL}, coverage profiles will be estimated for all transcripts in the
reference catalog.}

\item{nCores}{Integer, number of cores to use for parallel computations.}

\item{tx2gene}{A \code{data.frame} with transcript-to-gene mapping. Must have
at least two columns, named \code{tx} and \code{gene}.}

\item{genome}{A \code{BSgenome} object.}

\item{verbose}{Logical, whether to print progress messages.}
}
\value{
A list with one element per transcript. Each element is a list with
  five components:
\describe{
\item{\code{predCov}:}{The predicted coverage profile for the transcript.}
\item{\code{strand}:}{The annotated strand for the transcript.}
\item{\code{junctionCov}:}{A \code{data.frame} with genomic coordinates and
predicted read coverages for each junction in the transcript.}
\item{\code{aveFragLength}:}{The estimated average fragment length.}
\item{\code{note}:}{Either 'covOK', 'covError' or 'covNA'. If 'covOK',
\code{alpine} was able to predict a coverage profile. If 'covError' or
'covNA', the coverage profile could not be predicted, most likely because the
transcript is shorter than the fragment length or because no reads in the BAM
file overlapped the transcript. In both these cases, we impose a uniform
coverage profile for the transcript.}
}
}
\description{
This function predicts coverage profiles for each transcript in a specified
set of genes (or all transcripts in the reference catalog), based on the
fragment bias model estimated by \code{alpine}.
}
\examples{
\dontrun{
gtf <- system.file("extdata/Homo_sapiens.GRCh38.90.chr22.gtf.gz",
                   package = "jcc")
bam <- system.file("extdata/reads.chr22.bam", package = "jcc")
biasMod <- fitAlpineBiasModel(gtf = gtf, bam = bam,
                              organism = "Homo_sapiens",
                              genome = Hsapiens, genomeVersion = "GRCh38",
                              version = 90, minLength = 230,
                              maxLength = 7000, minCount = 10,
                              maxCount = 10000, subsample = TRUE,
                              nbrSubsample = 30, seed = 1, minSize = NULL,
                              maxSize = 220, verbose = TRUE)
tx2gene <- readRDS(system.file("extdata/tx2gene.sub.rds", package = "jcc"))
predCovProfiles <- predictTxCoverage(biasModel = biasMod$biasModel,
                                     exonsByTx = biasMod$exonsByTx,
                                     bam = bam, tx2gene = tx2gene,
                                     genome = Hsapiens,
                                     genes = c("ENSG00000070371",
                                               "ENSG00000093010"),
                                     nCores = 1, verbose = TRUE)
}

}
\references{
Soneson C, Love MI, Patro R, Hussain S, Malhotra D, Robinson MD: A junction
coverage compatibility score to quantify the reliability of transcript
abundance estimates and annotation catalogs. bioRxiv doi:10.1101/378539
(2018)

Love MI, Hogenesch JB, Irizarry RA: Modeling of RNA-seq fragment sequence
bias reduces systematic errors in transcript abundance estimation. Nature
Biotechnology 34(12):1287-1291 (2016).
}
\author{
Charlotte Soneson, Michael I Love
}
