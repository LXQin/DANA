#' Mean-Standard Deviation Plot
#'
#' Plots the marker-specific standard deviation over marker-specific means
#' for log2-normalized read counts. Bounds for the definition of positive
#' and negative controls can be highlighted by vertical lines (if given).
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#' @param tZero \emph{Optional.} Lower bound for negative controls.
#' @param tPoor \emph{Optional.} Upper bound for negative controls.
#' @param tWell \emph{Optional.} Lower bound for positive controls.
#' @param title Plot title
#' @param xlim.max Upper limit of the x-axis
#' @param ylim.max Upper limit of the y-axis
#'
#' @return ggplot object of the mean-sd plot
#' @export plotMeanSD
#'
#' @examples
#' counts <- matrix(rnbinom(1000, mu=3, size=0.01), nrow=50)
#' plotMeanSD(counts)
plotMeanSD <- function(raw, tZero, tPoor, tWell, title, xlim.max=NA, ylim.max=NA) {
  df <- data.frame(data.mean= rowMeans(log2(raw + 1)),
                   data.var = apply(log2(raw + 1), 1, stats::sd))

  if(missing(tZero)) {
    vline.t.zero <- ggplot2::geom_blank()
  } else {
    vline.t.zero <- ggplot2::geom_vline(xintercept = log2(tZero+1), color="#5851b8", linetype="twodash")
  }
  if(missing(tPoor)) {
    vline.t.poor <- ggplot2::geom_blank()
  } else {
    vline.t.poor <- ggplot2::geom_vline(xintercept = log2(tPoor+1), color="#5851b8", linetype="twodash")
  }
  if(missing(tWell)) {
    vline.t.well <- ggplot2::geom_blank()
  } else {
    vline.t.well <- ggplot2::geom_vline(xintercept = log2(tWell+1), color="#E7298A",  linetype="longdash")
  }

  if(missing(title)) {
    title <- ggplot2::geom_blank()
  } else {
    title <- ggplot2::ggtitle(title)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x=data.mean,y=data.var)) +
    ggplot2::geom_point(alpha=.25) +
    ggplot2::xlab("Mean") +
    ggplot2::ylab("Standard deviation") +
    vline.t.zero +
    vline.t.poor +
    vline.t.well +
    title +
    ggplot2::xlim(0, xlim.max) +
    ggplot2::ylim(0, ylim.max) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme_classic()
  return(p)
}



#' Read count histogram
#'
#' Plots a histogram of the marker-specific mean of log2-transformed
#' read counts. Bounds for the definition of positive
#' and negative controls can be highlighted by vertical lines (if given).
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#' @param tZero \emph{Optional.} Lower bound for negative controls.
#' @param tPoor \emph{Optional.} Upper bound for negative controls.
#' @param tWell \emph{Optional.} Lower bound for positive controls.
#' @param title Plot title
#' @param xlim.max Upper limit of the x-axis
#' @param ylim.max Upper limit of the y-axis
#'
#' @return ggplot object of the read count histogram
#' @export plotCountHist
#'
#' @examples
#' counts <- matrix(rnbinom(1000, mu=3, size=0.01), nrow=50)
#' plotCountHist(counts)
plotCountHist <- function(raw, binwidth=0.1, tZero, tPoor, tWell, title) {
  if(missing(tZero)) {
    vline.t.zero <- ggplot2::geom_blank()
  } else {
    vline.t.zero <- ggplot2::geom_vline(xintercept = log2(tZero+1), color="#5851b8", linetype="twodash")
  }
  if(missing(tPoor)) {
    vline.t.poor <- ggplot2::geom_blank()
  } else {
    vline.t.poor <- ggplot2::geom_vline(xintercept = log2(tPoor+1), color="#5851b8", linetype="twodash")
  }
  if(missing(tWell)) {
    vline.t.well <- ggplot2::geom_blank()
  } else {
    vline.t.well <- ggplot2::geom_vline(xintercept = log2(tWell+1), color="#E7298A",  linetype="longdash")
  }

  if(missing(title)) {
    title <- ggplot2::geom_blank()
  } else {
    title <- ggplot2::ggtitle(title)
  }

  df <- data.frame(log.expression = log2(rowMeans(raw) + 1))
  p <- ggplot2::ggplot(df, ggplot2::aes(x=log.expression)) +
    ggplot2::geom_histogram(binwidth = binwidth, color="black", fill="black") +
    vline.t.zero +
    vline.t.poor +
    vline.t.well +
    title +
    ggplot2::theme_classic()
  return(p)
}



#' Scatter plot DANA assessment metrics
#'
#' Plot function for a scatter plot of the DANA result metrics generated
#' using the DANA normalization assessment \code{\link{assessNormalization}}.
#' The DANA Assessment metrics for given normalized counts (for each
#' normalized count matrix) are provided as a \code{data.frame} with
#' the columns:
#' \describe{
#'   \item{cc}{\code{cc} measures the preservation of biological signals
#'   before versus after normalization.
#'   A high value indicates a high preservation of biological signals
#'   (\code{cc} <= 1)}
#'   \item{mscr}{\code{mscr} measures the relative reduction of handling before
#'   versus after normalization. A high \code{mscr} indicates higher
#'   removal of handling effects.}
#' }
#' When selecting a normalization method for the \code{raw} data, one should
#' aim for the best possible trade-off of hight cc and high mscr.
#' This function generates a \code{ggplot} scatter plot of the DANA metrics
#' for the selection of a most suitable normalization method.
#'
#' @param metrics \code{data.frame} of DANA metrics containing
#' the metrics \code{cc} and \code{mscr} as columns.
#' Each row represents a normalization method.
#' @param label.size \emph{Optional.} Text size of the plot labels.
#' @param label.repel \emph{Optional.} Logical if the
#' \code{ggrepel} package is used to space point labels.
#' Requires the package \code{ggrepel} to be installed.
#' @param fixed.limits.x \emph{Optional.} Logical if the limits of the x-axis
#' (mscr-) are fixed to [0,1]
#' @param fixed.limits.y \emph{Optional.} Logical if the limits of the y-axis
#' (cc+) are fixed to [0,1]
#'
#'
#' @return \code{ggplot} object of a scatter plot for the given DANA metrics.
#'
#' @export plotDANA
plotDANA <- function(metrics, label.size=3, label.repel=FALSE, fixed.limits.x=TRUE, fixed.limits.y=TRUE) {
  if (label.repel) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop("Package \"ggrepel\" needed for this function to work. Please install it.",
           call. = FALSE
      )
    }
  }

  scipen.prev <- getOption("scipen")
  options(scipen=4) # force non-scientific notation of x axis

  df <- data.frame(method=rownames(metrics),
                   cc = metrics$cc,
                   mscr = metrics$mscr)

  p <- ggplot2::ggplot(df, ggplot2::aes(x=mscr, y=cc, label=method)) +
       ggplot2::geom_point(alpha=1)

  if(label.repel) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = method),
                                      size=label.size, max.overlaps = Inf)
  } else {
    p <- p + ggplot2::geom_text(size=label.size, hjust=0, vjust=0)
  }

  p <- p +
    ggplot2::theme_classic() +
    ggplot2::xlab("mscr-; Relative reduction of handling effects") +
    ggplot2::ylab("cc+; Biological signal preservation")

  if(fixed.limits.x) {
    p <- p + ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0,1))
  } else {
    p <- p + ggplot2::scale_x_continuous(labels = scales::percent)
  }

  if(fixed.limits.y) {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1))
  } else {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent)
  }
  options(scipen=scipen.prev)  # reset scipen option
  return(p)
}



#' Correlation plot for partial correlations within polycistronic clusters
#'
#' This function generates a correlation (heatmap-type) plot for
#' partial correlations. Polycistronic clusters of the markers are
#' highlighted in the lower tri-diagonal part of the matrix while
#' partial correlations are shown in the upper tri-diagonal part.
#' Partial correlations can be computed using the \code{\link{partialCor}}
#' function.
#'
#' @param pCor Matrix of partial correlations. Rows and columns must be
#'   named by their corresponding markers
#' @param clusters Named Vector of clusters. Associates each marker
#'   to a polycistronic cluster. Usually generated using
#'   the function \code{\link{defineClusters}}.
#' @param title \emph{Optional.} Plot title
#'
#' @return \code{ggplot} object of a correlation plot for the given
#' partial correlations.
#' @export plotCor
#'
plotCor <- function(pCor, clusters, title){
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package \"ggnewscale\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }

  # plot title
  if(missing(title)) {
    title <- ggplot2::geom_blank()
  } else {
    title <- ggplot2::ggtitle(title)
  }

  # Convert gene-wise clusters to list of clusters
  geneClusters <- clusters
  clusters <- listClusters(clusters)

  # Set corr below diagonal to NA
  pCor[lower.tri(pCor)] <- NA
  # Set diagonal to NA
  diag(pCor) <- NA

  # Correlation data for plot
  data <- expand.grid(X=rownames(pCor), Y=colnames(pCor))[upper.tri(pCor, diag=FALSE),]
  data$Z <- as.vector(pCor[upper.tri(pCor, diag=FALSE)])

  clustered <- pCor
  clustered[] <- 0 # Set all to false

  diag(clustered) <- 1
  # color the miRNA clusters below the diagonal
  for (clust in clusters) {
    # find genes in clust in corr
    idx <- na.omit(match(unlist(clust), colnames(clustered)))
    # set edges in clusters to brigth color in the corrplot
    clustered[idx, idx][lower.tri(clustered[idx,idx])] <- 1
  }

  # Clustering data for plot
  data.clust <- expand.grid(X=rownames(clustered), Y=colnames(clustered))[lower.tri(clustered, diag=TRUE),]
  data.clust$Z <- as.vector(clustered[lower.tri(clustered, diag=TRUE)])

  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data=data.clust, ggplot2::aes(X, Y, fill=factor(Z))) +
    ggplot2::scale_fill_manual(name="Clusters",
                               values = c("0"="white", "1"="black"),
                               labels = c("un-clustered", "clustered")) +
    title +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(data=data, ggplot2::aes(X, Y, fill= Z)) +
    ggplot2::scale_x_discrete(limits = rev(levels(as.factor(data$X)))) +
    # ggplot2::scale_fill_gradient2(name = "Cor", midpoint=0,
    #                               low = "blue", mid = "white",
    #                               high = "red", limit = c(-1, 1),
    #                               space = "Lab") +
    ggplot2::scale_fill_gradientn(colours = c("blue4", "blue", "white", "red", " red4"),
                                  name = "Cor", limit = c(-1, 1),
                                  space = "Lab") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank()) +
    ggplot2::theme(aspect.ratio=1)
  return(p)
}





