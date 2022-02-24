#' Definition of negative and positive control miRNAs for DANA normalization assessment
#'
#'
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#'   The matrix rows and columns must be named, with rownames=miRNAs
#'   and colnames=samples.
#' @param tZero Lower (read count) bound for the definition of negative controls.
#' @param tPoor Upper (read count) bound for the definition of negative controls.
#' @param tWell Lower (read count) bound for the definition of positive controls.
#' @param clusters Named Vector of clusters. Associates each miRNA
#'   in \code{raw} to a polycistronic cluster. Usually generated using
#'   the function \code{\link{defineClusters}}.
#'
#' @return List of two vectors:
#' \describe{
#'   \item{posControls}{Vector of positive control markers (miRNAs) in \code{raw}.}
#'   \item{negControls}{Vector of negative control markers (miRNAs) in \code{raw}.}
#' }
#'
#' @export
#'
defineControls <- function(raw,
                           tZero,
                           tPoor,
                           tWell,
                           clusters) {

  ## Check Input
  # data matrices must have rownames
  if (is.null(rownames(raw))) {
    stop("Matrix raw must have rownames corresponding to miRNA names.")
  }


  genes <- rownames(raw)

  # convert gene-wise clusters to nested list
  geneClusters <- as.factor(clusters)
  clusters <- listClusters(geneClusters)


  ## Define negative control miRNAs
  negControls <- c()
  data.means <- rowMeans(raw)
  for(clust in clusters) {
    clust.genes <- unlist(clust)
    clust.genes.range <- names(which((data.means[clust.genes] <= tPoor) &
                                       (data.means[clust.genes] >= tZero)))
    if(length(clust.genes.range)!=0) {
      negControls <- c(negControls, clust.genes.range[length(clust.genes.range)])
    }
  }
  # cat("Number of negative control markers with RC in [", tZero, ", ",
  #     tPoor, "]: ", length(negControls), "\n" , sep="")

  ## Define positive control miRNAs
  posControls <- c()
  data.means <- rowMeans(raw)
  for(clust in clusters) {
    clust.genes <- unlist(clust)
    clust.genes.range <- names(which(data.means[clust.genes] >= tWell))
    if(length(clust.genes.range) >= 2) {
      posControls <- c(posControls, clust.genes.range)
    }
  }
  # cat("Number of positive control markers with RC in [", tWell, ", inf): ",
  #     length(posControls), "\n" , sep="")

  # # miRNAs in clusters with 2 or more members
  # names(clusters)[which(clusters %in% names(table(clusters))[table(clusters) >= 2])]

  return(list(
    negControls = negControls,
    posControls = posControls
  ))
}



#' Polycistronic clustering of miRNAs
#'
#' Defines polycistronic clusters for a given set of miRNAs based on their
#' location on the genome.
#'
#' @param genes List of genes for which the miRNA clustering is computed.
#' @param chr Vector that specifies the chromosome for each miRNA in
#'   \code{genes}.
#' @param pos Vector of the base-pair position of each miRNA in \code{genes} on
#'   its chromosome, which is specified in \code{chr}. This position could be
#'   e.g. the start, end or averaged ((start+end)/2) base-pair position on the
#'   associated chromosome.
#' @param threshold Maximum distance threshold in base-pairs for the difintion
#'   of polycistronic clusters.
#'
#' @return Vector of clusters. Associates each miRNA in \code{genes} to a
#'   polycistronic cluster.
#' @export
defineClusters <- function(genes, chr, pos, threshold=10000) {
  if ((length(genes) != length(chr)) || length(genes) != length(pos)){
    stop("Lengths of vectors genes, chr, and pos must agree.")
  }
  # Check for NAs
  if(anyNA(genes)) stop("genes array contains NAs.")
  if(anyNA(chr)) stop("chr array contains NAs.")
  if(anyNA(pos)) stop("pos array contains NAs.")

  names(chr) <- names(pos) <- genes
  chr <- as.factor(chr)
  clusters <- rep(NA, length(genes))
  names(clusters) <- genes

  # loop over the chromosomes and order/sort the genes on each chromosome by their location
  for (curChr in levels(chr)) {
    # find genes on chromosome curChr
    genesOnChrm <- genes[curChr == chr]
    # locally sort genes on chromosome curChr
    localOrdering <- order(pos[curChr == chr])
    genesOnChrm <- genesOnChrm[localOrdering]

    # add first gene to first cluster
    curCluster <- genesOnChrm[1]
    clusters[genesOnChrm[1]] <- curCluster

    # sequentially add more genes
    for (i in 2:length(genesOnChrm)) {
      if (abs(pos[genesOnChrm[i]] - pos[genesOnChrm[i-1]]) <= threshold) {
        # add gene to current cluster
        clusters[genesOnChrm[i]] <- curCluster
      } else {
        # new cluster
        curCluster <- genesOnChrm[i]
        clusters[genesOnChrm[i]] <- curCluster
      }
    }
  }
  return(clusters)
}
