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

sdrs <- function(x, grouping, dims, prec.est = NULL, ...){

  classes <- unique(grouping)
  data <- lapply(1:length(classes), function(i){
    x[grouping == classes[[i]], ]
  })

  n <- lapply(data, nrow)
  S <- lapply(data, stats::cov)
  xbar <- lapply(data, colMeans)

  if(is.list(prec.est)) S_inv <- prec.est

  if(is.null(prec.est)) {
    S_inv <- lapply(S, function(x) solve(x, tol = NULL))
  } else if(prec.est == "Haff"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      Haffshrinkage(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "RidgeShrinkage"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      wang_est(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "BGP16"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      p <- ncol(S[[i]])
      HDShOP::InvCovShrinkBGP16(n[[i]], p, diag(p), solve(S[[i]], tol = NULL))$S
    })
  } else if(prec.est == "SCPME"){
    S_inv <- SCPME_qda(x, grouping, ...)
    lam <- sapply(seq_along(S_inv), function(i) S_inv[[i]]$Tuning[2])
    S_inv <- lapply(seq_along(S_inv), function(i) S_inv[[i]]$Omega)
  } else if(prec.est == "glasso"){
    # glasso_wrapper <- function(S,n,  lam = NULL, ...){
    #   out <- vector(mode = "list", length = length(S))
    #   if(is.null(lam)) lam <- rep(NULL, length(S))
    #   for (i in seq_along(S)) {
    #     out[[i]] <- glasso::glasso(s = S[[i]],
    #                                rho = lam[[i]], nobs = n[[i]])$wi
    #   }
    #   out
    # }
    # S_inv <- glasso_wrapper(data, ...)

    glasso_out <- lapply(data, function(x){
      CVglasso::CVglasso(x, trace = "none", crit.cv = "loglik")
    })
    S_inv <- lapply(glasso_out, function(x) x$Omega)
    lam <- sapply(glasso_out, function(x) x$Tuning[[2]])
  }



  projectedMeanDiffs <- do.call(cbind,
                                lapply(2:length(xbar), function(i) {
                                  S_inv[[i]] %*% xbar[[i]] -
                                    S_inv[[1]] %*% xbar[[1]]
                                }))

  Sdiffs <- do.call(cbind,
                    lapply(2:length(S), function(i) {
                      S[[i]] - S[[1]]
                    }))

  M <- cbind(projectedMeanDiffs, Sdiffs)
  U <- as.matrix(svd(M)$u)
  #Return Projection Matrix
  out <- list("M" = M, "ProjectionMatrix" = t(U[,dims]))

  if(!is.null(prec.est)){
    if(prec.est == "SCPME" | prec.est == "glasso") out$lam <-  lam
  }
  out
}

sdrs2 <- function(x, grouping, dims, prec.est = "glasso", ...){

  classes <- unique(grouping)
  data <- lapply(1:length(classes), function(i){
    x[grouping == classes[[i]], ]
  })

  n <- lapply(data, nrow)
  S <- lapply(data, stats::cov)
  xbar <- lapply(data, colMeans)
  xbarbar <- colMeans(x)

  if(is.list(prec.est)) S_inv <- prec.est

  if(is.null(prec.est)) {
    S_inv <- lapply(S, function(x) solve(x, tol = NULL))
  } else if(prec.est == "Haff"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      Haffshrinkage(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "RidgeShrinkage"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      wang_est(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "BGP16"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      p <- ncol(S[[i]])
      HDShOP::InvCovShrinkBGP16(n[[i]], p, diag(p), solve(S[[i]], tol = NULL))$S
    })
  } else if(prec.est == "SCPME"){
    S_inv <- SCPME_qda(x, grouping, ...)
    lam <- sapply(seq_along(S_inv), function(i) S_inv[[i]]$Tuning[2])
    S_inv <- lapply(seq_along(S_inv), function(i) S_inv[[i]]$Omega)
  } else if(prec.est == "glasso"){
    # glasso_wrapper <- function(S,n,  lam = NULL, ...){
    #   out <- vector(mode = "list", length = length(S))
    #   if(is.null(lam)) lam <- rep(NULL, length(S))
    #   for (i in seq_along(S)) {
    #     out[[i]] <- glasso::glasso(s = S[[i]],
    #                                rho = lam[[i]], nobs = n[[i]])$wi
    #   }
    #   out
    # }
    # S_inv <- glasso_wrapper(data, ...)

    # glasso_out <- lapply(data, function(x){
    #   CVglasso::CVglasso(x, trace = "none", crit.cv = "loglik")
    # })



    glasso_wrapper <- function(data, S, lam = NULL){
      out <- vector(mode = "list", length = length(S))
      for (i in seq_along(S)) {
        S_biased <- S[[i]]*(n[[i]] - 1)/n[[i]]
        temp <- CVglasso::CVglasso(data[[i]], S = S_biased,
                                              trace = "none", crit.cv = "loglik", lam = lam[[i]])
        if(is.null(lam)){
          min.ratio <- 0.01
          while(temp$Lambdas[1] == temp$Tuning[[2]] && min.ratio >= 0.000625){
            min.ratio <- min.ratio/2
            temp <- CVglasso::CVglasso(x, trace = "none", crit.cv = "loglik", lam.min.ratio = min.ratio)
          }
        }
        out[[i]] <- temp
      }
      out
    }

    glasso_out <- glasso_wrapper(data, S, ...)
    S <- lapply(glasso_out, function(x) x$Sigma)
    S_inv <- lapply(glasso_out, function(x) x$Omega)
    lam <- sapply(glasso_out, function(x) x$Tuning[[2]])

    # glasso_out <- lapply(seq_along(data), function(i){
    #   S_biased <- S[[i]]*(n[[i]] - 1)/n[[i]]
    #   # lams <- eigen(S_biased)$values
    #   temp <- CVglasso::CVglasso(data[[i]], S = S_biased,
    #                              trace = "none", crit.cv = "loglik"
    #                              # ,
    #                              # lam = sort(lams)
    #                              )

    # })

  }

  projectedMeanDiffs <- do.call(cbind,
                                lapply(1:length(xbar), function(i) {
                                  S_inv[[i]] %*% (xbar[[i]] - xbarbar)
                                    # S_inv[[1]] %*% xbar[[1]]
                                }))

  Sdiffs <- do.call(cbind,
                    lapply(2:length(S), function(i) {
                      S[[i]] - S[[1]]
                    }))

  Sdiff_svd <- svd(Sdiffs)

  top <- lik(Sdiff_svd$d^2)$max

  Sdiffs <- Sdiff_svd$u[,1:top]

  M <- cbind(normalize(projectedMeanDiffs), Sdiffs)
  # MMT <- M %*% t(M)
  # PM_cv <- PMA::PMD.cv(MMT, type = "ordered", sumabsus = seq(1, floor(sqrt(nrow(M))), by = 1), sumabss = NULL)
  # U <- PMA::PMD(MMT, type = "ordered", sumabsu = PM_cv$bestsumabsu, K = nrow(M), v = PM_cv$v.init, center = FALSE)$v
  # U <- PMA::SPC(M %*% t(M), K = nrow(M))$v
  U <- as.matrix(svd(M, nu = nrow(M))$u)
  # U <- eigen(M %*% t(M))$vectors
  # U_r <- as.matrix(svd(M, nu = nrow(M))$u)[,1:ncol(M)]
  # U <- U_r %*% t(U_r)

  #Return Projection Matrix
  out <- list("M" = M, "ProjectionMatrix" = t(U[,dims]))

  if(!is.null(prec.est)){
    if(prec.est == "SCPME" | prec.est == "glasso") out$lam <-  lam
  }
  out
}

lik <- function(d){
  p <- length(d)
  out <- vector(length = p)
  for(q in 1:(p-1)){
    d_list <- list(d[1:q], d[(q+1):p])
    mu <- lapply(d_list, mean)
    s2 <- lapply(d_list, function(x) ifelse(length(x) == 1, as.numeric(0), var(x)) )
    pooled_sd <- sqrt(((q-1)*s2[[1]] + (p-q-1)*s2[[2]] )/(p-2))

    out[[q]] <- sapply(seq_along(d_list), function(i){
      sum(dnorm(d_list[[i]], mu[[i]], pooled_sd, log = TRUE))
    }) |>
      sum()
  }
  out[[p]] <- sum(dnorm(d, mean(d), sqrt(var(d)), log = TRUE))
  ## Return
  list("lik" = out, "max" = which.max(out))
}

normalize <- function(x) {x / sqrt(sum(x^2))}
# dot <- function(x, y) {
#   stopifnot(length(x) == length(y))
#   sum(Conj(y) * x)
# }
#
# norm <- function(x) sqrt(dot(x, x))
#
# is.zero <- function(x) norm(x) < sqrt(.Machine$double.eps)
# normalize <- function(x) {
#   if (is.zero(x)) {
#     stop("Input `x` is numerically equivalent to the zero vector.", call. = FALSE)
#   } else {
#     x / norm(x)
#   }
# }

# SYS2 <- function(x, grouping, dims, ...){
#   classes <- unique(grouping)
#   data <- lapply(1:length(classes), function(i){
#     x[grouping == classes[[i]], ]
#   })
#
#   n <- lapply(data, nrow)
#   S <- lapply(data, stats::cov)
#   xbar <- lapply(data, colMeans)
#
#   # St_inv <- lapply(seq_along(S), function(i, ...){
#   #   Haffshrinkage(S = S[[i]], n = n[[i]], tol = NULL)
#   # })
#
#   St_inv <- lapply(seq_along(S), function(i, ...){
#     RidgeShrinkage(S = S[[i]], n = n[[i]], tol = NULL)
#   })
#
#   projectedMeanDiffs <- do.call(cbind,
#                                 lapply(2:length(xbar), function(i) {
#                                   St_inv[[i]] %*% xbar[[i]] -
#                                     St_inv[[1]] %*% xbar[[1]]
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
#   t(U[,dims])
# }

