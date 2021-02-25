#' Partial correlation estimation
#'
#' Partial correlation estimation using  the neighborhood selection method (see
#' Meinshausen, Buehlmann (2006)) for the estimation of the inverse covariance
#' matrix (precision matrix)
#'
#' @references
#' \href{https://projecteuclid.org/euclid.aos/1152540754}{Neighborhood
#' Selection}
#'
#' @param X Data matrix. Size (samples x parameters): n x p. Columns should be
#' named corresponding to marker/miRNA names.
#' @param scale Logical if the data is scaled column-wise to mean 0 and var 1.
#'
#' @return Estimated partial correlation matrix. Shape p x p
#' @export
#'
partialCor <- function(X, scale=TRUE) {
  ## Pre-processing

  X <- t(X)
  n <- dim(X)[1]
  p <- dim(X)[2]

  # remove genes with read count 0
  X <- X[, !(colSums(abs(X))==0)]
  # Remove genes with constant read count across all markers
  X <- X[, apply(X, 2, function(y) length(unique(y)) > 1)]
  cat("Partial correlation estimation. Filtered",
      p-dim(X)[2], "genes from the data\n")
  p <- dim(X)[2]

  # modified log2 transform
  X <- log2(X+1)

  # scale the data by centering and dividing by standard deviation
  if(scale) {
    X <- apply(
      X, 2, function (y) (y-mean(y)) / stats::sd(y) ^ as.logical(stats::sd(y))
    )
  }


  ## Estimate precision using the neighborhood selection method
  ## The tuning parameter is selected using Bayesian Information Criterion

  betas <- matrix(-1, p, p)  # augmented betas
  tuning <- rep(NA, p)  # tuning parameters
  # compute betas using the lasso estimator with Bayesian information criterion
  for (i in 1:p) {
    lasso.fit <- glmnet::glmnet(x=X[, -i], y=X[, i], intercept=FALSE)
    # BIC criterion
    tLL <- lasso.fit$nulldev - stats::deviance(lasso.fit)  # log-likelihood
    dof <- lasso.fit$df  # degrees-of-freedom
    lasso.BIC <- log(n)*dof - tLL
    BIC.optimal <- which.min(lasso.BIC)
    betas[-i, i] <- lasso.fit$beta[, BIC.optimal]
    tuning[i] <- lasso.fit$lambda[BIC.optimal]
  }

  # precision matrix = inverse covriance matrix
  Theta <- matrix(0, p, p)

  # compute diagonal elements
  for (i in 1:p) {
    Theta[i, i] <- n / ( sum((X[, i] - X[, -i] %*% betas[-i, i]) ^ 2) )
  }

  # compute off-diagonal elements
  for (i in 1:p) {
    for (j in 1:p) {
      Theta[i, j] <- -0.5 * (Theta[i, i] * betas[j, i] +
                               Theta[j, j] * betas[i, j])
    }
  }

  colnames(Theta) <- rownames(Theta) <- colnames(X)
  partial.correlations <- prec2corr(Theta)
}



#' Pearson correlation estimation
#'
#' Estimaetes pearson correlations for a given data matrix.
#'
#'
#' @param X Data matrix. Size (samples x parameters): n x p
#'
#' @return Estimated Pearson correlation matrix. Shape p x p
pearsonCor <- function(X) {
  ## Pre-processing

  X <- t(X)
  n <- dim(X)[1]
  p <- dim(X)[2]

  # remove genes with read count 0
  X <- X[, !(colSums(abs(X))==0)]
  # Remove genes with constant read count across all markers
  X <- X[, apply(X, 2, function(y) length(unique(y)) > 1)]
  cat("Pearson correlation estimation. Filtered",
      p-dim(X)[2], "genes from the data\n")
  p <- dim(X)[2]

  # modified log2 transform
  X <- log2(X+1)
  pearson.correlations <- stats::cor(X, method="pearson")
}




#' Computes partial correlation matrix from a given precision matrix
#' @param P A square precision matrix
#' @return Partial correlation matrix computed from P
prec2corr <- function(P) {

  # Test if matrix is square
  if(dim(P)[1] != dim(P)[2]) {
    stop("Error: Function prec2cor expected a square matrix.")
  }

  n <- dim(P)[1]  # dimension of matrix
  C <-matrix(1, nrow=n, ncol=n)

  # compute the partial correlations
  for (i in 1:n) {
    for (j in 1:n) {
      C[i,j] <- -P[i,j] / sqrt(P[i,i] * P[j,j])
    }
  }

  # catch any correlations >1 or <-1 and cap the values
  C[C > 1] <- 1
  C[C < (-1)] <- -1

  rownames(C) <- rownames(P)
  colnames(C) <- colnames(P)

  return(C)
}
