#' Polycistronic clustering of miRNAs
#'
#' Defines polycistronic clusters for a given set of miRNAs based on their
#' chromosomal location. For microRNAs on the human genome (homo sapiens) using
#' mirBAse v22 notation, chromosomal coordinates are provided. These can be used
#' e.g. for data from TCGA. If you use miRNA annotation different from miRBase
#' v22, you must provide your own chromosomal location including chromosome
#' ('chr') and nucleotide position on the chromosome ('pos').
#'
#' @param genes List of genes for which the miRNA clustering is computed.
#' @param chr \emph{Optional.} Vector that specifies the chromosome for each
#'   miRNA in \code{genes}. If not given, miRBase v22 chromosomal coordinates is
#'   used.
#' @param pos \emph{Optional.} Vector of the nucleotice (base-pair) position of
#'   each miRNA in \code{genes} on its chromosome, which is specified in
#'   \code{chr}. This position could be e.g. the start, end or averaged
#'   ((start+end)/2) base-pair position on the associated chromosome. If not
#'   given, miRBase v22 chromosomal coordinates is used.
#' @param threshold Maximum distance threshold in base-pairs for the difintion
#'   of polycistronic clusters. Default is 10,000, which corresponds to
#'   miRBase's cluster definition.
#'
#' @return Vector of clusters. Associates each miRNA in \code{genes} to a
#'   polycistronic cluster.
#' @export
defineClusters <- function(genes, chr, pos, threshold=10000,...) {
  # Check for NAs
  if(anyNA(genes)) stop("genes array contains NAs.")
  # Check if only one of chr and pos is missing
  if((missing(chr) + missing(pos)) == 1){
    stop("You must provide either both chr and pos, or neither. If neither are given, miRBase coordinates will be used.")
  }
  # Check if mirbase coordinates should be used (if chr and pos are missing)
  if(missing(chr) && missing(pos)) {
    use.mirbase <- TRUE
  } else {
    use.mirbase <- FALSE
  }

  # use all lower-case miRNA names
  genes <- tolower(genes)

  # Set up chr and pos
  if(use.mirbase) {
    # check for genes that are not present in the coordinate data frame
    genes.in.coords <- rownames(miRNA.coordinates)[na.omit(match(genes, rownames(miRNA.coordinates)))]
    if(length(genes)-length(genes.in.coords) > 1) {
      cat(length(genes)-length(genes.in.coords),
          " either have no unique chromosomal coordinate or could not be found in the hsa miRBase coordinate set. Reducing gene set to ",
          length(genes.in.coords), "genes.\n")
    }
    if(length(genes.in.coords) == 0) {
      stop("No miRNAs could be found in the hsa miRBase coordinate set. Please verify that your data set uses miRBase miRNA notation.")
    }
    genes <- genes.in.coords

    # set chr and pos from miRBase coordinates
    chr <- DANA::miRNA.coordinates[genes, "chromosome"]
    pos <- DANA::miRNA.coordinates[genes, "position"]
  } else {  # use provided chr and pos for clustering
    # Check for same length of genes and chr and pos
    if ((length(genes) != length(chr)) || length(genes) != length(pos)){
      stop("Lengths of vectors genes, chr, and pos must agree.")
    }
  }
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
    if(length(genesOnChrm) > 1) {
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
  }
  return(clusters)
}




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
  # Check ranges of tZero, tPoor, tWell
  if(tZero < 0) {
    stop("tZero must be positive. Negative control genes are defined using the intervall [tZero, tPoor].")
  }
  if(tZero > tPoor) {
    stop("tPoor must be greater than tZero. Negative control genes are defined using the intervall [tZero, tPoor].")
  }
  if(tWell < tPoor) {
    stop("tWell must be greater than tPoor.")
  }

  genes <- rownames(raw)



  # convert gene-wise clusters to nested list
  geneClusters <- as.factor(clusters)
  clusters <- listClusters(geneClusters)


  ## Define negative control miRNAs
  negControls <- c()
  data.means <- rowMeans(raw)
  # Clustered genes
  for(clust in clusters) {
    clust.genes <- unlist(clust)
    clust.genes.range <- names(which((data.means[clust.genes] <= tPoor) &
                                       (data.means[clust.genes] >= tZero)))
    if(length(clust.genes.range)!=0) {
      # Add a single gene from cluster clust
      negControls <- c(negControls, clust.genes.range[length(clust.genes.range)])
    }
  }
  # Un-clustered genes
  unclustered.genes <- genes[!(genes %in% names(geneClusters))]
  genes.range <- names(which((data.means[unclustered.genes] <= tPoor) &
                               (data.means[unclustered.genes] >= tZero)))
  if(length(genes.range)!=0) {
    negControls <- c(negControls, genes.range)
  }

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


