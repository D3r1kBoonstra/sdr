#### Haff Shrinkage Estimator ####
Haffshrinkage <- function(X, S = NULL, n = NULL, ...){
  ## T(U_i)
  tu <- function(S, n, p) {
    U <- p * det(S) ^ (1 / p) / sum(diag(S))
    min((4 * (p ^ 2 - 1) / ((n - p - 2) * p ^ 2)), 1) * U ^ (1 / p)
  }
  # tu <- function(S, n, p) {
  #   U <- p * sum(eigen(S, only.values = TRUE)$values) ^ (1 / p) / sum(diag(S))
  #   min((4 * (p ^ 2 - 1) / ((n - p - 2) * p ^ 2)), 1) * U ^ (1 / p)
  # }

  ## If X is given
  if(is.null(S)){
    ### X as matrix
    if(!is.matrix(X)) X <- is.matrix(X)
    ### Cov of X
    S <- stats::cov(X)
    ### nrow of X
    n <- nrow(X)
  }
  ## If S is given
  else {
    ### S as matrix
    if(!is.matrix(S)) S <- as.matrix(S)
    ### Checking if n is given
    if(is.null(n)) stop("If using S, then n must be given")
  }
  p <- ncol(S)
  S_inv <- solve(S, ...)
  tu <- tu(S, n, p)
  ## Computing Haff Shrinkage estimator
  (1 - tu)*(n - p - 2)*S_inv +
    ((tu * (n*p - p - 2)) / sum(diag(S)))*diag(1, p)
}
