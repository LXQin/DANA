% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{norm.PoissonSeq}
\alias{norm.PoissonSeq}
\title{PoissonSeq Normalization}
\usage{
norm.PoissonSeq(raw)
}
\arguments{
\item{raw}{Raw read count matrix (rows = genes, cols = samples).}

\item{groups}{Sample groups (subtypes)}
}
\value{
A list with the following elements:
\describe{
  \item{dataNormalized}{Normalized read counts.}
  \item{scalingFactor}{Normalizing factors for scaling the raw data.}
}
}
\description{
PoissonSeq normalization for miRNA-Seq data.
}
\examples{
rawCounts <- matrix(0:20, nrow = 7)
normCounts <- norm.PoissonSeq(rawCounts)
}
