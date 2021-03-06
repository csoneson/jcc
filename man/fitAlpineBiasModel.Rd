% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitAlpineBiasModel.R
\name{fitAlpineBiasModel}
\alias{fitAlpineBiasModel}
\title{Fit fragment bias model with alpine}
\usage{
fitAlpineBiasModel(gtf, bam, organism, genome, genomeVersion, version,
  minLength = 600, maxLength = 7000, minCount = 500, maxCount = 10000,
  subsample = TRUE, nbrSubsample = 200, seed = 1, minSize = NULL,
  maxSize = NULL, verbose = FALSE)
}
\arguments{
\item{gtf}{Path to gtf file with genomic features. Preferably in Ensembl
format.}

\item{bam}{Path to bam file with read alignments to the genome.}

\item{organism}{The organism (e.g., 'Homo_sapiens'). This argument will be
passed to \code{ensembldb::ensDbFromGtf}.}

\item{genome}{A \code{BSgenome} object.}

\item{genomeVersion}{Genome version (e.g., 'GRCh38'). This argument will be
passed to \code{ensembldb::ensDbFromGtf}.}

\item{version}{The version of the reference annotation (e.g., 90). This
argument will be passed to \code{ensembldb::ensDbFromGtf}.}

\item{minLength, maxLength}{Minimum and maximum length of single-isoform genes
used to fit fragment bias model.}

\item{minCount, maxCount}{Minimum and maximum read coverage of single-isoform
genes used to fit fragment bias model.}

\item{subsample}{Whether to subsample the set of single-isoform genes
satisfying the \code{minLength}, \code{maxLength}, \code{minCount} and
\code{maxCount} criteria before fitting the fragment bias model.}

\item{nbrSubsample}{If \code{subsample} is \code{TRUE}, the number of genes
to subsample.}

\item{seed}{If \code{subsample} is \code{TRUE}, the random seed to use to
ensure reproducibility.}

\item{minSize, maxSize}{Smallest and largest fragment size to consider. One or
both of these can be \code{NULL}, in which case it is estimated as the 2.5
or 97.5 percentile, respectively, of estimated fragment sizes in the
provided data.}

\item{verbose}{Logical, whether to print progress messages.}
}
\value{
A list with three elements: \describe{ \item{\code{biasModel}:}{The
  fitted fragment bias model.} \item{\code{exonsByTx}:}{A \code{GRangesList}
  object with exons grouped by transcript.} \item{\code{transcripts}:}{A
  \code{GRanges} object with all the reference transcripts.} }
}
\description{
This function provides a wrapper around some of the functions from the
\code{alpine} package. Given a gtf file and a bam file with reads aligned to
the genome, it will find single-isoform genes (with lengths and expression
levels within given ranges) and use the observed read coverages to fit a
fragment bias model.
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
}

}
\references{
Soneson C, Love MI, Patro R, Hussain S, Malhotra D, Robinson MD:
A junction coverage compatibility score to quantify the reliability of
transcript abundance estimates and annotation catalogs. bioRxiv
doi:10.1101/378539 (2018).

Love MI, Hogenesch JB, Irizarry RA: Modeling of RNA-seq fragment sequence
bias reduces systematic errors in transcript abundance estimation. Nature
Biotechnology 34(12):1287-1291 (2016).
}
\author{
Charlotte Soneson, Michael I Love
}
