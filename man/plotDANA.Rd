% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotDANA}
\alias{plotDANA}
\title{Scatter plot DANA assessment metrics}
\usage{
plotDANA(metrics, label.size = 3, label.repel = FALSE)
}
\arguments{
\item{metrics}{\code{data.frame} of DANA metrics containing
the metrics \code{cc} and \code{mcr} as columns.
Each row represents a normalization method.}

\item{label.size}{\emph{Optional.} Text size of the plot labels.}

\item{label.repel}{\emph{Optional.} Logical if the
\code{ggrepel} package is used to space point labels.
Requires the package \code{ggrepel} to be installed.}
}
\value{
\code{ggplot} object of a scatter plot for the given DANA metrics.
}
\description{
Plot function for a scatter plot of the DANA result metrics generated
using the DANA normalization assessment \code{\link{assessNormalization}}.
The DANA Assessment metrics for given normalized counts (for each
normalized count matrix) are provided as a \code{data.frame} with
the columns:
\describe{
  \item{cc}{\code{cc} measures the preservation of biological signals
  before versus after normalization.
  A high value indicates a high preservation of biological signals
  (\code{cc} <= 1)}
  \item{mcr}{\code{mcr} measures the relative reduction of handling before
  versus after normalization. A high \code{mcr} indicates higher
  removal of handling effects.}
}
When selecting a normalization method for the \code{raw} data, one should
aim for the best possible trade-off of hight cc and high mcr.
This function generates a \code{ggplot} scatter plot of the DANA metrics
for the selection of a most suitable normalization method.
}
