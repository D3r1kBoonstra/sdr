library("SCPME")

# qda_shrink <- function(x, grouping, prec.est = NULL, xtest = NULL, ...) {
#   x <- as.matrix(x)
#   classes <- unique(grouping)
#   data <- lapply(1:length(classes), function(i) {
#     x[grouping == classes[[i]],]
#   })
#   p <- ncol(x)
#   n <- lapply(data, nrow)
#   N <- length(grouping)
#   priors <- lapply(n, function(x)
#     x / N)
#   xbar <- lapply(data, colMeans)
#   S <- lapply(data, stats::cov)
#   if(is.null(prec.est)) {
#     S_inv <- lapply(S, function(x) solve(x, tol = NULL))
#   } else if(prec.est == "Haff"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       Haffshrinkage(S = S[[i]], n = n[[i]], tol = NULL)
#     })
#   } else if(prec.est == "RidgeShrinkage"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       RidgeShrinkage(S = S[[i]], n = n[[i]], tol = NULL)
#     })
#   } else if(prec.est == "BGP16"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       p <- ncol(S[[i]])
#       HDShOP::InvCovShrinkBGP16(n[[i]], p, diag(p), solve(S[[i]], tol = NULL))$S
#     })
#   } else if(prec.est == "SCPME"){
#     S_inv <- SCPME_qda(x, grouping, ...)
#     S_inv <- lapply(seq_along(S_inv), function(i) S_inv[[i]]$Omega)
#   }
#
#   if (!is.null(xtest)) {
#     x <- as.matrix(xtest)
#     N <- nrow(x)
#   }
#
#   d <- lapply(seq_along(data), function(k) {
#     # logs <- log(det(S[[i]])) - 2*log(priors[[i]])
#     out <- vector(length = N)
#     for (i in 1:N) {
#       diff <- as.matrix((x[i, ] - xbar[[k]]))
#       out[[i]] <- -.5 * (t(diff) %*% S_inv[[k]] %*% diff - determinant(S_inv[[k]])$mod) + log(priors[[k]])
#     }
#     out
#   }) |>
#     do.call(cbind, args = _)
#
#   classes[apply(d, 1, which.max)]
# }

SCPME_qda <- function(x, grouping, gamma = NULL, lam = NULL, type = "L1", standardize_xbar = FALSE, ...){
  if(is.null(gamma)) gamma <- rep(1, length(unique(grouping)))
  x <- as.matrix(x)
  classes <- unique(grouping)
  data <- lapply(1:length(classes), function(i) {
    x[grouping == classes[[i]],]
  })
  p <- ncol(x)
  xbar <- lapply(data, colMeans)
  S <- lapply(data, cov)
  out <- vector("list", length(data))

  for (k in seq_along(data)) {
    if(type == "L1"){
      A <-  diag(p)
      B <-  diag(p)
      C <-  matrix(0, nrow = p, ncol = p)
    } else if (type == "qda"){
      standardize <- function(x) {(x - mean(x))/sd(x)}
      if(standardize_xbar) xbar <- lapply(xbar, standardize)
      B <- cbind(matrix(xbar[[k]], ncol = 1), gamma[[k]]*diag(p))
      A <- t(B)
      C <- matrix(0, nrow = nrow(A), ncol = ncol(B))
    }
    if(is.null(lam)) lam <- rep(NULL, length(data))
    out[[k]] <-  SCPME::shrink(X = data[[k]],
                               A = A, B = B, C = C,
                               crit.cv = "loglik", trace = "none", lam = lam[[k]], ...)
  }
  out
}

