% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{defineControls}
\alias{defineControls}
\title{Definition of negative and positive control miRNAs for DANA normalization assessment}
\usage{
defineControls(raw, tZero, tPoor, tWell, clusters)
}
\arguments{
\item{raw}{Raw read count matrix (rows = genes, cols = samples).
The matrix rows and columns must be named, with rownames=miRNAs
and colnames=samples.}

\item{tZero}{Lower (read count) bound for the definition of negative controls.}

\item{tPoor}{Upper (read count) bound for the definition of negative controls.}

\item{tWell}{Lower (read count) bound for the definition of positive controls.}

\item{clusters}{Named Vector of clusters. Associates each miRNA
in \code{raw} to a polycistronic cluster. Usually generated using
the function \code{\link{defineClusters}}.}
}
\value{
List of two vectors:
\describe{
  \item{posControls}{Vector of positive control markers (miRNAs) in \code{raw}.}
  \item{negControls}{Vector of negative control markers (miRNAs) in \code{raw}.}
}
}
\description{
Definition of negative and positive control miRNAs for DANA normalization assessment
}
