#' Applies normalization methods to given data set
#'
#' Applies single or multiple normalization methods to the given data.
#'   Available methods are:
#'   \itemize{
#'     \item Total Count (TC)
#'     \item Upper Quartile (UQ)
#'     \item Median (median)
#'     \item Trimmed Median of Means (TMM)
#'     \item DESeq
#'     \item PoissonSeq
#'     \item Quantile Normalization (QN)
#'     \item Remove Unwanted Variation (RUVg, RUVr, and RUVs)
#'   }
#'
#' @param data Un-normalized count data set of shape n x p, where n is the
#'   number of samples and p is the number of markers.
#' @param groups Vector of length n that maps the samples to sample groups.
#' @param method Vector of normalization methods that are applied to the data.
#'   Available methods are: \code{c("TC", "UQ", "median", "TMM", "DESeq",
#'   "PoissonSeq", "QN", "RUV")}. Select one or multiple methods. By default all
#'   normalization methods will be applied.
#'
#' @export
#' @return List of normalized read counts for selected normalization methods.
applyNormalization <- function(data,
                               groups,
                               method = c(
                                 "TC", "UQ", "Median", "TMM", "DESeq",
                                 "PoissonSeq", "QN", "RUV"
                               )) {
  method <- match.arg(
    method,
    c("TC", "UQ", "median", "TMM", "DESeq", "PoissonSeq", "QN", "RUV"),
    several.ok = TRUE
  )

  data.normalized <- list()

  if ("TC" %in% method) {
    data.normalized <- append(data.normalized, list("TC" = norm.TC(data, groups)$dataNormalized))
  }
  if ("UQ" %in% method) {
    data.normalized <- append(data.normalized, list("UQ" = norm.UQ(data, groups)$dataNormalized))
  }
  if ("median" %in% method) {
    data.normalized <- append(data.normalized, list("Med" = norm.Med(data, groups)$dataNormalized))
  }
  if ("TMM" %in% method) {
    data.normalized <- append(data.normalized, list("TMM" = norm.TMM(data, groups)$dataNormalized))
  }
  if ("DESeq" %in% method) {
    data.normalized <- append(data.normalized, list("DESeq" = norm.DESeq(data, groups)$dataNormalized))
  }
  if ("PoissonSeq" %in% method) {
    data.normalized <- append(data.normalized, list("PoissonSeq" = norm.PoissonSeq(data)$dataNormalized))
  }
  if ("QN" %in% method) {
    data.normalized <- append(data.normalized, list("QN" = norm.QN(data)$dataNormalized))
  }
  if ("RUV" %in% method) {
    data.normalized <- append(data.normalized, list("RUVg" = norm.RUV(data, groups, method = "RUVg")$dataNormalized))
    data.normalized <- append(data.normalized, list("RUVs" = norm.RUV(data, groups, method = "RUVs")$dataNormalized))
    data.normalized <- append(data.normalized, list("RUVr" = norm.RUV(data, groups, method = "RUVr")$dataNormalized))
  }

  # if ("SVA" %in% method) {
  #   data.normalized <- append(data.normalized, "SVA" = list(norm.SVA(data, groups)$dataNormalized))
  # }

  return(data.normalized)
}



# https://academic.oup.com/bib/article/14/6/671/189645
#' Total Count (TC) Normalization
#'
#' Total count normalization for miRNA-Seq data.
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#' @param groups Sample groups (subtypes)
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{dataNormalized}{Normalized read counts.}
#'   \item{scalingFactor}{Normalizing factors for scaling the raw data.}
#' }
#' @export
#'
#' @examples
#' rawCounts <- matrix(0:20, nrow = 7)
#' normCounts <- norm.TC(rawCounts)
norm.TC <- function(raw, groups = rep(1, ncol(raw))) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  dat.DGE <- edgeR::DGEList(
    counts = matrix(raw, ncol = length(groups)),
    group = factor(groups),
    genes = rownames(raw)
  )
  scalingFactor <- dat.DGE$samples$lib.size / 1e6
  dataNormalized <- t(t(raw) / scalingFactor)
  # scalingFactor <- dat.DGE$samples$lib.size / mean(dat.DGE$samples$lib.size)
  # dataNormalized <- round(t(t(raw)/scalingFactor))
  return(list(
    dataNormalized = dataNormalized,
    scalingFactor = scalingFactor
  ))
}



# calcNormFactors could also obtain the UQ normalization factor
# We compute cor() between its UQ and our UQ normalized benchmark data
# and find it to be 1
#' Upper Quartile (UQ) Normalization
#'
#' Upper quartile (UQ) normalization for miRNA-Seq data.
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#' @param groups Sample groups (subtypes)
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{dataNormalized}{Normalized read counts.}
#'   \item{scalingFactor}{Normalizing factors for scaling the raw data.}
#' }
#' @export
#'
#' @examples
#' rawCounts <- matrix(0:20, nrow = 7)
#' normCounts <- norm.UQ(rawCounts)
norm.UQ <- function(raw, groups = rep(1, ncol(raw))) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  dat.DGE <- edgeR::DGEList(
    counts = matrix(raw, ncol = length(groups)),
    group = factor(groups),
    genes = rownames(raw)
  )
  q.factor <- apply(dat.DGE$counts, 2, function(x) stats::quantile(x[x != 0], probs = 0.75))
  scalingFactor <- q.factor / 1e6
  dataNormalized <- t(t(raw) / scalingFactor)
  # scalingFactor <- q.factor/mean(dat.DGE$samples$lib.size)
  # dataNormalized <- round(t(t(raw)/scalingFactor))
  return(list(
    dataNormalized = dataNormalized,
    scalingFactor = scalingFactor
  ))
}



# https://academic.oup.com/bib/article/14/6/671/189645
#' Median (Med) Normalization
#'
#' Median (Med) normalization for miRNA-Seq data.
#'
#' @param raw Raw read count matrix (rows = genes, cols = samples).
#' @param groups Sample groups (subtypes)
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{dataNormalized}{Normalized read counts.}
#'   \item{scalingFactor}{Normalizing factors for scaling the raw data.}
#' }
#' @export
#'
#' @examples
#' rawCounts <- matrix(0:20, nrow = 7)
#' normCounts <- norm.Med(rawCounts)
norm.Med <- function(raw, groups = rep(1, ncol(raw))) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  dat.DGE <- edgeR::DGEList(
    counts = matrix(raw, ncol = length(groups)),
    group = factor(groups),
    genes = rownames(raw)
  )
  m.factor <- apply(dat.DGE$counts, 2, function(x) stats::median(x[x != 0]))
  scalingFactor <- m.factor / 1e6
  dataNormalized <- t(t(raw) / scalingFactor)
  # scalingFactor <- m.factor/mean(dat.DGE$samples$lib.size)
  # dataNormalized <- round(t(t(raw)/scalingFactor))
  return(list(
    dataNormalized = dataNormalized,
    scalingFactor = scalingFactor
  ))
}


## Trimmed Mean of M-values (TMM)

# https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# page 14
# cpm() could compute normalized count per million based on DGE result
# We compute cor() between cpm() and our normalized benchmark data and find
# it to be 1
norm.TMM <- function(raw, groups = rep(1, ncol(raw))) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  dat.DGE <- edgeR::DGEList(
    counts = matrix(raw, ncol = length(groups)),
    group = factor(groups),
    genes = rownames(raw)
  )
  d <- edgeR::calcNormFactors(dat.DGE, method = "TMM")
  scalingFactor <- d$samples$norm.factors * d$samples$lib.size / 1e6
  dataNormalized <- t(t(raw) / scalingFactor)
  # scalingFactor <- d$samples$norm.factors * d$samples$lib.size /
  #   mean(dat.DGE$samples$lib.size)
  # dataNormalized <- round(t(t(raw)/scalingFactor))
  return(list(
    dataNormalized = dataNormalized,
    scalingFactor = scalingFactor,
    dge.normed = d
  ))
}





## DESeq

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf
# page 4
norm.DESeq <- function(raw, groups = rep(1, ncol(raw))) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package \"DESeq\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  if (!requireNamespace("BiocGenerics", quietly = TRUE)) {
    stop("Package \"BiocGenerics\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  condition <- data.frame(SampleName = colnames(raw), Condition = factor(groups))
  rownames(condition) = colnames(raw)
  dat.DGE <- DESeq2::DESeqDataSetFromMatrix(countData = raw, colData = condition, design = ~ Condition)
  dat.DGE <- DESeq2::estimateSizeFactors(dat.DGE)
  scalingFactor <- DESeq2::sizeFactors(dat.DGE)
  dataNormalized <- BiocGenerics::counts(dat.DGE, normalized = T)
  return(list(dataNormalized = dataNormalized,
              scalingFactor = scalingFactor))
}


## PoissonSeq

# https://cran.r-project.org/web/packages/PoissonSeq/PoissonSeq.pdf
norm.PoissonSeq <- function(raw) {
  if (!requireNamespace("PoissonSeq", quietly = TRUE)) {
    stop("Package \"PoissonSeq\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  scalingFactor <- PoissonSeq::PS.Est.Depth(raw)
  dataNormalized <- t(t(raw) / scalingFactor)
  return(list(
    dataNormalized = dataNormalized,
    scalingFactor = scalingFactor
  ))
}


## Quantile Normalization (QN)

# http://jtleek.com/genstats/inst/doc/02_05_normalization.html
norm.QN <- function(raw, filter = FALSE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package \"limma\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  if (filter == TRUE) {
    raw <- log2(raw + 1)
    raw <- raw[rowMeans(raw) > 2, ]
  } else {
    raw <- log2(raw + 1)
  }
  dat.log.normed <- limma::normalizeQuantiles(as.matrix(raw))
  dataNormalized <- 2^dat.log.normed - 1
  colnames(dataNormalized) <- colnames(raw)
  rownames(dataNormalized) <- rownames(raw)
  return(list(dataNormalized = dataNormalized))
}


#' Normalization By Remove Unwanted Variation Using Replicate Samples (RUVs)
#' Normalize the dataset using RUV using replicate samples, and return the normalized dataset with adjusting factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#' @return list, containing \code{dat.normed} (normalized dataset), and the \code{adjust.factor} (adjusting factors) for the design matrix. The normalized dataset could only used for exploration, and adjusting factors are recommended as a covariate in the downstream analysis.
#'
#' @export
#'
#' @references \href{http://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf}{RUVSeq Tutorial}
#'
#' @examples
#' test.RUVs <- norm.RUVs(data.test, data.group)
norm.RUV <- function(raw, groups, method = c("RUVg", "RUVs", "RUVr")) {
  if (!require("Biobase")) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("RUVSeq", quietly = TRUE)) {
    stop("Package \"RUVSeq\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("EDASeq", quietly = TRUE)) {
    stop("Package \"EDASeq\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }

  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  genes <- rownames(dat.ruv)
  condition <- factor(groups)
  set <- EDASeq::newSeqExpressionSet(
    as.matrix(dat.ruv),
    phenoData = data.frame(condition, row.names = colnames(dat.ruv))
  )
  design <- stats::model.matrix(
    ~condition,
    data = data.frame(condition,
      row.names = colnames(dat.ruv)
    )
  )
  y <- edgeR::DGEList(counts = dat.ruv, group = condition)
  y <- edgeR::calcNormFactors(y, method = "upperquartile")
  y <- edgeR::estimateGLMCommonDisp(y, design)
  y <- edgeR::estimateGLMTagwiseDisp(y, design)
  fit <- edgeR::glmFit(y, design)
  lrt <- edgeR::glmLRT(fit, coef = 2)
  top <- edgeR::topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:0.15 * nrow(raw)]))]

  if (method == "RUVg") {
    t <- RUVSeq::RUVg(set, spikes, k = 1)
    dataNormalized <- EDASeq::normCounts(t)
    return(list(
      dataNormalized = dataNormalized,
      adjust.factor = t$W
    ))
  } else if (method == "RUVs") {
    differences <- RUVSeq::makeGroups(condition)
    controls <- rownames(dat.ruv)
    t <- RUVSeq::RUVs(set, controls, k = 1, differences)
    dataNormalized <- EDASeq::normCounts(t)
    return(list(
      dataNormalized = dataNormalized,
      adjust.factor = t$W
    ))
  } else if (method == "RUVr") {
    design <- stats::model.matrix(~condition, data = Biobase::pData(set))
    y <- edgeR::DGEList(counts = dat.ruv, group = condition)
    y <- edgeR::calcNormFactors(y, method = "upperquartile")
    y <- edgeR::estimateGLMCommonDisp(y, design)
    y <- edgeR::estimateGLMTagwiseDisp(y, design)
    fit <- edgeR::glmFit(y, design)
    res <- stats::residuals(fit, type = "deviance")
    setUQ <- EDASeq::betweenLaneNormalization(set, which = "upper")
    controls <- rownames(dat.ruv)
    t <- RUVSeq::RUVr(setUQ, controls, k = 1, res)
    dataNormalized <- EDASeq::normCounts(t)

    return(list(
      dataNormalized = dataNormalized,
      adjust.factor = t$W
    ))
  }
}



## SVA for Sequencing (SVASeq)

# https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
# page 13
# calibrator genes include tag "cali"
# Normalized data only used for exploration
# https://www.biostars.org/p/121489/
norm.SVA <- function(raw, groups) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package \"sva\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.sva <- raw[filter, ]
  genes <- rownames(dat.sva)
  mod1 <- stats::model.matrix(~groups)
  mod0 <- cbind(mod1[, 1])
  dat0 <- as.matrix(dat.sva)
  svseq <- sva::svaseq(dat0, mod1, mod0, n.sv = 1)$sv
  adjust <- cbind(mod1, svseq)
  hat <- solve(t(adjust) %*% adjust) %*% t(adjust)
  beta <- (hat %*% t(raw))
  P <- ncol(mod1)
  dataNormalized <- raw - t(as.matrix(adjust[, -c(1:P)]) %*% beta[-c(1:P), ])
  return(list(
    dataNormalized = dataNormalized,
    adjust.factor = svseq
  ))
}
