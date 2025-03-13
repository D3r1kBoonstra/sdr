wang_est <- function(X, S = NULL, n = NULL, ...){

  if(is.null(S)){
    ### X as matrix
    if(!is.matrix(X)) X <- is.matrix(X)
    ### Cov of X
    S <- cov(X)
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

  L <- function(lambda){
    out <- vector(length = length(lambda))
    for (i in 1:length(lambda)) {
      tr <- sum(diag((1/lambda[[i]])*S + diag(p)))
      a1 <-  1 - (1/p)*tr^(-1)
      a2 <- (1/p)*tr^(-1) - (1/p)*tr^(-2)
      yhat <- p/n
      R1 <- a1/(1 - yhat*a1)
      R2 <- a1/((1 - yhat*a1)^(3)) - a2/((1 - yhat*a1)^(4))

      out[[i]] <- 1 - R1^2/R2
    }
    # Return
    out
  }
  alph <- function(Beta){
    tr <- sum(diag((1/Beta)*S + diag(p)))
    a1 <-  1 - (1/p)*tr^(-1)
    a2 <- (1/p)*tr^(-1) - (1/p)*tr^(-2)
    yhat <- p/n
    R1 <- a1/(1 - yhat*a1)
    R2 <- a1/((1 - yhat*a1)^(3)) - a2/((1 - yhat*a1)^(4))
    # Return
    R1/R2
  }
  eigs <- eigen(cov(S))$values
  Beta <- optimize(L, c(min(eigs), max(eigs)))$minimum
  alpha <- alph(Beta)

  # Return
  alpha*solve(S + Beta*diag(p))
}

haff_est <- function(X, S = NULL, n = NULL, ...){
  tu <- function(S, n, p) {
    U <- p * det(S) ^ (1 / p) / sum(diag(S))
    min((4 * (p ^ 2 - 1) / ((n - p - 2) * p ^ 2)), 1) * U ^ (1 / p)
  }

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

# sdrs <- function(x, grouping, dims, prec.est = "glasso", ...){
#
#   classes <- unique(grouping)
#   data <- lapply(1:length(classes), function(i){
#     x[grouping == classes[[i]], ]
#   })
#
#   n <- sapply(data, nrow)
#   n_order <- order(n, decreasing = TRUE)
#   n <- n[n_order]
#   prior <- lapply(seq_along(n), function(i) n[[i]]/sum(n))
#   data <- data[n_order]
#   S <- lapply(data, stats::cov)
#   xbar <- lapply(data, colMeans)
#
#   # if(is.list(prec.est)) S_inv <- prec.est
#
#   if(is.null(prec.est)) {
#     S_inv <- lapply(S, function(x) solve(x, tol = NULL))
#   } else if(prec.est == "Haff"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       haff_est(S = S[[i]], n = n[[i]], tol = NULL)
#     })
#   } else if(prec.est == "Wang"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       wang_est(S = S[[i]], n = n[[i]], tol = NULL)
#     })
#   } else if(prec.est == "Bodnar"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       p <- ncol(S[[i]])
#       HDShOP::InvCovShrinkBGP16(n[[i]], p, diag(p), solve(S[[i]], tol = NULL))$S
#     })
#   } else if(prec.est == "MRY"){
#     S_inv <- SCPME_qda(x, grouping, ...)
#     lam <- sapply(seq_along(S_inv), function(i) S_inv[[i]]$Tuning[2])
#     S_inv <- lapply(seq_along(S_inv), function(i) S_inv[[i]]$Omega)
#   } else if(prec.est == "glasso"){
#     glasso_wrapper <- function(data, S, lam = NULL){
#       out <- vector(mode = "list", length = length(S))
#
#       S_biased <- lapply(seq_along(S), function(i) S[[i]]*(n[[i]] - 1)/n[[i]])
#       # if(is.null(lam)){
#       #   lam <- lapply(S_biased, function(x){
#       #     eig_vals <- eigen(x)$values
#       #     max_off_diag <- max(abs(x[upper.tri(x)]))
#       #     exp(seq(
#       #       log(min(eig_vals)),
#       #       log(min(max(eig_vals), max_off_diag)),
#       #       length.out = 10)
#       #     )
#       #   })
#       # }
#
#       for (i in seq_along(S)) {
#
#         out[[i]] <- CVglasso::CVglasso(data[[i]], S = S_biased[[i]],
#                                        trace = "none", crit.cv = "loglik", lam = lam[[i]])
#         # if(is.null(lam)){
#         #   min.ratio <- 0.01
#         #   while(temp$Lambdas[1] == temp$Tuning[[2]] && min.ratio >= 0.000625){
#         #     min.ratio <- min.ratio/2
#         #     temp <- CVglasso::CVglasso(x, trace = "none", crit.cv = "loglik", lam.min.ratio = min.ratio)
#         #   }
#         # }
#         # out[[i]] <- temp
#       }
#       out
#     }
#     glasso_out <- glasso_wrapper(data, S, ...)
#     # S <- lapply(glasso_out, function(x) x$Sigma)
#     S_inv <- lapply(glasso_out, function(x) x$Omega)
#     lam <- sapply(glasso_out, function(x) x$Tuning[[2]])
#   }
#
#   projectedMeanDiffs <- do.call(cbind,
#                                 lapply(2:length(xbar), function(i) {
#                                   S_inv[[i]] %*% xbar[[i]] -
#                                     S_inv[[1]] %*% xbar[[1]]
#                                 }))
#
#   Sdiffs <- do.call(cbind,
#                     lapply(2:length(S), function(i) {
#                       S[[i]] - S[[1]]
#                     }))
#
#   M <- cbind(projectedMeanDiffs, Sdiffs)
#   U <- as.matrix(svd(M)$u)
#   #Return Projection Matrix
#   out <- list("M" = M, "ProjectionMatrix" = t(U[,dims]))
#
#   if(!is.null(prec.est)){
#     if(prec.est == "SCPME" | prec.est == "glasso") out$lam <-  lam
#   }
#   out
# }
#
# sdrs2 <- function(x, grouping, dims, prec.est = "glasso", ...){
#
#   classes <- unique(grouping)
#   data <- lapply(1:length(classes), function(i){
#     x[grouping == classes[[i]], ]
#   })
#
#   n <- sapply(data, nrow)
#   n_order <- order(n, decreasing = TRUE)
#   n <- n[n_order]
#   prior <- lapply(seq_along(n), function(i) n[[i]]/sum(n))
#   data <- data[n_order]
#   S <- lapply(data, stats::cov)
#   # S_common <- cov(x)
#   xbar <- lapply(data, colMeans)
#   xbarbar <- colMeans(x)
#
#   # if(is.list(prec.est)) S_inv <- prec.est
#
#   if(is.null(prec.est)) {
#     S_inv <- lapply(S, function(x) solve(x, tol = NULL))
#   } else if(prec.est == "Haff"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       Haffshrinkage(S = S[[i]], n = n[[i]], tol = NULL)
#     })
#   } else if(prec.est == "Wang"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       wang_est(S = S[[i]], n = n[[i]], tol = NULL)
#     })
#   } else if(prec.est == "Bodnar"){
#     S_inv <- lapply(seq_along(S), function(i, ...){
#       p <- ncol(S[[i]])
#       HDShOP::InvCovShrinkBGP16(n[[i]], p, diag(p), solve(S[[i]], tol = NULL))$S
#     })
#   } else if(prec.est == "MRY"){
#     S_inv <- SCPME_qda(x, grouping, ...)
#     lam <- sapply(seq_along(S_inv), function(i) S_inv[[i]]$Tuning[2])
#     S_inv <- lapply(seq_along(S_inv), function(i) S_inv[[i]]$Omega)
#   } else if(prec.est == "glasso"){
#
#     glasso_wrapper <- function(data, S, lam = NULL){
#       out <- vector(mode = "list", length = length(S))
#
#       S_biased <- lapply(seq_along(S), function(i) S[[i]]*(n[[i]] - 1)/n[[i]])
#       # if(is.null(lam)){
#       #   lam <- lapply(S_biased, function(x){
#       #     eig_vals <- eigen(x)$values
#       #     eig_zero <- sapply(eig_vals, is.zero)
#       #     eig_vals <- eig_vals[!eig_zero]
#       #     max_off_diag <- max(abs(x[upper.tri(x)]))
#       #     exp(seq(
#       #       log(min(eig_vals)),
#       #       log(min(max(eig_vals), max_off_diag)),
#       #       length.out = 10)
#       #     )
#       #   })
#       # }
#
#       for (i in seq_along(S)) {
#
#         out[[i]] <- CVglasso::CVglasso(data[[i]], S = S_biased[[i]],
#                                               trace = "none", crit.cv = "loglik", lam = lam[[i]])
#         # if(is.null(lam)){
#         #   min.ratio <- 0.01
#         #   while(temp$Lambdas[1] == temp$Tuning[[2]] && min.ratio >= 0.000625){
#         #     min.ratio <- min.ratio/2
#         #     temp <- CVglasso::CVglasso(x, trace = "none", crit.cv = "loglik", lam.min.ratio = min.ratio)
#         #   }
#         # }
#         # out[[i]] <- temp
#       }
#       out
#     }
#     glasso_out <- glasso_wrapper(data, S, ...)
#     # S <- lapply(glasso_out, function(x) x$Sigma)
#     S_inv <- lapply(glasso_out, function(x) x$Omega)
#     lam <- sapply(glasso_out, function(x) x$Tuning[[2]])
#
#   }
#
#   # S <- lapply(S_inv, function(x) solve(x))
#
#   projectedMeanDiffs <- do.call(cbind,
#                                 lapply(1:length(xbar), function(i) {
#                                   S_inv[[i]] %*% (xbar[[i]] - xbarbar)
#                                 }))
#
#   Sdiffs <- do.call(cbind,
#                     lapply(2:length(S), function(i) {
#                       S[[i]] - S[[1]]
#                     }))
#
#   Sdiff_svd <- svd(Sdiffs)
#
#   top <- lik(Sdiff_svd$d^2)$max
#
#   # out <- core(Sigmas = S, ns = n, numdir.test = TRUE, numdir = ncol(S[[1]]))
#   # d <- which.min(out$bic)
#   # Sdiffs <- out$Gammahat[[d]]
#
#   # Sdiifs <- with(Sdiff_svd,
#   #                as.matrix(u[,1:top]) %*% as.matrix(diag(d))[1:top, 1:top])
#
#   Sdiffs <- Sdiff_svd$u[,1:top]
#
#   M <- cbind(apply(projectedMeanDiffs, 2, normalize), Sdiffs)
#   U <- as.matrix(svd(M, nu = nrow(M))$u)
#
#   #Return Projection Matrix
#   out <- list("M" = M, "ProjectionMatrix" = t(U[,dims]))
#
#   if(!is.null(prec.est)){
#     if(prec.est == "SCPME" | prec.est == "glasso") out$lam <-  lam
#   }
#   out
# }
#
# lik <- function(d){
#   p <- length(d)
#   out <- vector(length = p)
#   for(q in 1:(p-1)){
#     d_list <- list(d[1:q], d[(q+1):p])
#     mu <- lapply(d_list, mean)
#     s2 <- lapply(d_list, function(x) ifelse(length(x) == 1, as.numeric(0), var(x)) )
#     pooled_sd <- sqrt(((q-1)*s2[[1]] + (p-q-1)*s2[[2]] )/(p-2))
#
#     out[[q]] <- sapply(seq_along(d_list), function(i){
#       sum(dnorm(d_list[[i]], mu[[i]], pooled_sd, log = TRUE))
#     }) |>
#       sum()
#   }
#   out[[p]] <- sum(dnorm(d, mean(d), sqrt(var(d)), log = TRUE))
#   ## Return
#   list("lik" = out, "max" = which.max(out))
# }
#
# dot <- function(x, y) {
#   stopifnot(length(x) == length(y))
#   sum(Conj(y) * x)
# }
# norm <- function(x) sqrt(dot(x, x))
# is.zero <- function(x) norm(x) < sqrt(.Machine$double.eps)
#
# normalize <- function(x) {
#   # if (is.zero(x)) {
#   #   stop("Input `x` is numerically equivalent to the zero vector.", call. = FALSE)
#   # } else {
#     x / norm(x)
#   # }
# }
