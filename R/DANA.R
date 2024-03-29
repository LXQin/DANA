
#' Data-driven miRNA sequencing Normalization Assessment
#'
#'
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#'   The rows and columns of the count matrix must be named,
#'   where \code{rownames(raw)} are the marker names
#'   and \code{colnames(raw)} are the sample names.
#' @param normalized Named list of normalized count matrices.
#'   Each matrix holds the normalized read count matrix corresponding to a
#'   normalization method under study.
#'   Each list member must be named (e.g. after the used normalization).
#'   Each matrix in \code{normalized} must be named where
#'   the row names are the marker names
#'   and the column names are the sample names.
#'   A list of normalized counts can be generated using
#'   the \code{\link{applyNormalization}} function.
#' @param negControls Vector of negative control markers as generated by
#'   the function \code{\link{defineControls}}.
#' @param posControls Vector of positive control markers as generated by
#'   the function \code{\link{defineControls}}.
#' @param clusters Named Vector of clusters. Associates each miRNA
#'   in \code{raw} to a polycistronic cluster. Usually generated using
#'   the function \code{\link{defineClusters}}.
#'
#' @return DANA Assessment metrics for the provided normalized counts (for each
#'   normalized count matrix).
#'   DANA computes two assessment metrics:
#' \describe{
#'   \item{cc}{\code{cc} measures the preservation of biological signals
#'   before versus after normalization.
#'   A high value indicates a high preservation of biological signals
#'   (\code{cc} <= 1).
#'   In particular, \code{cc} is the concordance correlation coefficient of the
#'   within-cluster partial correlation among positive controls before and
#'   after normalization.}
#'   \item{mscr}{\code{mscr} measures the relative reduction of handling before
#'   versus after normalization. A high \code{mscr} indicates higher
#'   removal of handling effects.
#'   In particular, \code{mscr} is the mean-squared correlation reduction in
#'   negative controls before and after normalization.}
#' }
#' When selecting a normalization method for the \code{raw} data, one should
#' aim for the best possible trade-off of hight cc and high mscr.
#' @export
#'
assessNormalization <- function(raw,
                                normalized,
                                negControls,
                                posControls,
                                clusters) {
  ## Constants
  # Minimun number of controls (hard threshold)
  min.num.pos.controls <- 10
  min.num.neg.controls <- 10


  ## Check Input
  # data matrices must have rownames
  if (is.null(rownames(raw))) {
    stop("Matrix raw must have rownames corresponding to miRNA names.")
  }
  if (any(sapply(sapply(normalized, rownames), is.null))) {
    stop("All matrices in normalized must have rownames corresponding to miRNA names.")
  }
  # control miRNAs must be present in raw data
  if (!all(negControls %in% rownames(raw))) {
    stop("Could not find all negative control markers in raw data.")
  }
  if (!all(posControls %in% rownames(raw))) {
    stop("Could not find all positive control markers in raw data.")
  }
  # Markers in normalized must be present in raw
  for(normDat in normalized) {
    if(!all(rownames(normDat) %in% rownames(raw))) {
      stop("Normalized data contains markers that are not found in raw data.")
    }
  }
  # Clusters must give clustering information at least for all positive controls
  if (is.null(names(clusters))) {
    stop("Clusters must be named vector. Clustering information must contain
          at least all positive control clusters.")
  }
  if(!all(posControls %in% names(clusters))) {
    stop("Clustering information missing for some/all positive controls.
          Clustering information must contain be given at least all positive
          control clusters.")
  }


  ## Parameters
  # convert gene-wise clusters to nested list
  geneClusters <- as.factor(clusters)
  clusters <- listClusters(geneClusters)
  num.norm <- length(normalized)


  ### ESTIMATE CORRELATIONS ----------------------------------------------------


  ## Compute correlations in negative controls
  # raw data
  raw.corNeg <- pearsonCor(raw[negControls, ])
  # normalized data
  norm.corNeg <- list()
  for (i in 1:num.norm) {
    controlGenes <- negControls
    # Check if all negative control markers are present in normalized data set
    if (!all(negControls %in% rownames(normalized[[i]]))) {
      controlGenes <- intersect(rownames(normalized[[i]]), negControls)
      warning(
        paste0(length(negControls) - length(controlGenes),
               " negative control markers not found in normalized data set ",
               names(normalized)[i],
               ". Reducing number of negative controls for this data set from ",
               length(negControls), " to ", length(controlGenes))
      )
    }
    if (length(controlGenes) < min.num.neg.controls) {
      # Number of control markers too small
      stop(
        paste0("Number of negative controls (",
               length(controlGenes),
               ") too small for normalized data set ", names(normalized)[i])
      )
    }
    norm.corNeg <-
      append(norm.corNeg,
             list(pearsonCor(normalized[[i]][controlGenes, ])))
  }
  names(norm.corNeg) <- names(normalized)


  ## Compute correlations in positive controls
  # raw data
  raw.corPos <- partialCor(raw[posControls, ], scale=TRUE)
  # normalized data
  norm.corPos <- list()
  for (i in 1:num.norm) {
    controlGenes <- posControls
    # Check if all positive control markers are present in normalized data set
    if (!all(posControls %in% rownames(normalized[[i]]))) {
      controlGenes <- intersect(rownames(normalized[[i]]), posControls)
      warning(
        paste0(length(posControls) - length(controlGenes),
               " positive control markers not found in normalized data set",
               names(normalized)[i],
               ". Reducing number of positive controls for this data set to ",
               length(controlGenes))
      )
    }
    if (length(controlGenes) < min.num.pos.controls) {
      # Number of control markers too small
      stop(
        paste0("Number of positive control markers (",
               length(controlGenes),
               ") too small for normalized data set ", names(normalized)[i])
      )
    }
    norm.corPos <-
      append(norm.corPos,
             list(partialCor(normalized[[i]][controlGenes, ], scale=TRUE)))
  }
  names(norm.corPos) <- names(normalized)


  ### ASSESS NORMALIZATION -----------------------------------------------------

  ## Compute metric mscr- for negative controls
  mscr <- rep(NA, num.norm)
  for (i in 1:num.norm) {
    mscr[i] <- compute.mscr(raw.corNeg, norm.corNeg[[i]])
  }

  ## Compute metric cc+ for positive controls
  cc <- rep(NA, num.norm)
  for (i in 1:num.norm) {
    cc[i] <- compute.cc(raw.corPos, norm.corPos[[i]], clusters)
  }

  ## Return result metrics
  metrics <- data.frame(cc, mscr)
  rownames(metrics) <- names(normalized)
  return(metrics)
}



#' Compute the mean correlation reduction in negative controls
#' @keywords internal
compute.mscr <- function(rawCor, normCor) {
  varZero.rawCor <- sum(rawCor[upper.tri(rawCor)]^2) / sum(upper.tri(rawCor))
  varZero.normCor <- sum(normCor[upper.tri(normCor)]^2) / sum(upper.tri(normCor))
  return((varZero.rawCor - varZero.normCor) / varZero.rawCor)
}



#' Compute the concordance correlation of clustered positive controls
#' @keywords internal
compute.cc <- function(rawCor, normCor, clusters) {
  # remove genes that are not present in both models
  rawCor <- rawCor[!is.na(match(colnames(rawCor), colnames(normCor))),
                   !is.na(match(colnames(rawCor), colnames(normCor)))]
  normCor <- normCor[!is.na(match(colnames(normCor), colnames(rawCor))),
                     !is.na(match(colnames(normCor), colnames(rawCor)))]

  # consider only cluster with multiple genes
  clusters <- clusters[lengths(clusters) > 1]
  clusters.raw <- c()
  clusters.norm <- c()
  for (clust in clusters) {
    # only consider cluster genes that are positive controls
    clust.genes <- clust[stats::na.omit(match(colnames(rawCor), clust))]
    if (length(clust.genes) < 2) {
      next  # disregard clusters with less than 2 genes
    }
    # subset of correlations in the cluster
    rawCor.clust <- rawCor[clust.genes, clust.genes]
    clusters.raw <- c(clusters.raw, rawCor.clust[upper.tri(rawCor.clust)])
    normCor.clust <- normCor[clust.genes, clust.genes]
    clusters.norm <- c(clusters.norm, normCor.clust[upper.tri(normCor.clust)])
  }

  ## Compute the concordance correlation between nonzero partial correlations
  # idx <- clusters.raw | clusters.norm
  # if (sum(idx)>0) {
  #   cc <- DescTools::CCC(as.vector(clusters.raw[idx]),
  #                        as.vector(clusters.norm[idx]))$rho.c$est
  # } else {
  #   cc <- 0
  # }

  cc <- DescTools::CCC(as.vector(clusters.raw),as.vector(clusters.norm))$rho.c$est

  return(cc)
}












