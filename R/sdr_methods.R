sdr.fit.pca <- function(object, x, y, ytype, dims, ...){

  M <- S <- cov(x)
  N <- diag(ncol(M))


  eigs <- geigen::geigen(M, N)
  beta <- as.matrix(with(eigs, vectors[,ncol(vectors):1])[,dims])
  eigvalues <- with(eigs, values[length(values):1])[dims]
  ## Return
  out <- list("ProjectedData" = x %*% beta,
       "ProjectionMatrix" = beta,
       "eigvalues" = eigvalues,
       "M" = M, "N" = N,
       "x" = x, "y" = y,
       "ytype" = ytype
       )
  class(out) <- c("sdr")
  out
}

sdr.fit.sir <- function(object, x, y, ytype, dims, priors = NULL, ...){

  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }


  data_list <- data_list_fn(x, slices)
  if(is.null(priors)) priors <- priors_fn(data_list = data_list)

  xbarbar <- as.matrix(colMeans(x))
  N <- S <- cov(x)
  # S_inv_sqrt <- mat_power(S, -0.50)
  xbar <- lapply(data_list, function(x) as.matrix(colMeans(x)))


  M <- lapply(seq_along(xbar), function(i){
    diff <- xbar[[i]] - xbarbar
    priors[[i]] * diff %*% t(diff)
  }) |>
    Reduce(`+`, x = _)

  # beta <- S_inv_sqrt %*% eigen(S_inv_sqrt %*% M %*% S_inv_sqrt)$vectors
  eigs <- geigen::geigen(M, N)
  beta <- as.matrix(with(eigs, vectors[,ncol(vectors):1])[,dims])
  eigvalues <- with(eigs, values[length(values):1])[dims]
  ## Return
  out <- list("ProjectedData" = x %*% beta,
              "ProjectionMatrix" = beta,
              "eigvalues" = eigvalues,
              "M" = M, "N" = N,
              "x" = x, "y" = y,
              "ytype" = ytype, "slices" = slices
  )
  class(out) <- c("sdr")
  out
}

sdr.fit.save <- function(object, x, y, ytype, dims, priors = NULL, ...){
  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }
  p <- ncol(x)
  N <- S_x <- cov(x)
  S_inv_sqrt <- mat_power(S_x, -.5)
  S_sqrt <- mat_power(S_x, .5)
  z <- scale(x, center = TRUE, scale = FALSE) %*% S_inv_sqrt

  data_list <- data_list_fn(z, slices)
  if(is.null(priors)) priors <- priors_fn(data_list = data_list)

  S_zy <- lapply(data_list, cov)
  E <- lapply(seq_along(S_zy), function(i){
    priors[[i]]* mat_power(diag(p) - S_zy[[i]], 2)
  }) |>
    Reduce(`+`, x = _)

  M <- S_sqrt %*% E %*% S_sqrt

  # beta <- S_inv_sqrt %*% eigen(S_inv_sqrt %*% M %*% S_inv_sqrt)$vectors
  eigs <- geigen::geigen(M, N)
  beta <- as.matrix(with(eigs, vectors[,ncol(vectors):1])[,dims])
  eigvalues <- with(eigs, values[length(values):1])[dims]
  ## Return
  out <- list("ProjectedData" = x %*% beta,
              "ProjectionMatrix" = beta,
              "eigvalues" = eigvalues,
              "M" = M, "N" = N,
              "x" = x, "y" = y,
              "ytype" = ytype, "slices" = slices
  )
  class(out) <- c("sdr")
  out
}

sdr.fit.sir2 <- function(object, x, y, ytype, dims, priors = NULL, regularize = FALSE, lambda = NULL, ...){
  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }

  p <- ncol(x)
  N <- S_x <- cov(x)
  S_inv_sqrt <- mat_power(S_x, -.5)
  z <- scale(x, center = TRUE, scale = FALSE) %*% S_inv_sqrt

  data_list <- data_list_fn(z, slices)
  if(is.null(priors)) priors <- priors_fn(data_list = data_list)

  xbar <- lapply(data_list, function(x) as.matrix(colMeans(x)))
  M_sir <- lapply(seq_along(xbar), function(i){
    priors[[i]] * xbar[[i]] %*% t(xbar[[i]])
  }) |>
    Reduce(`+`, x = _)

  S <- lapply(data_list, cov)
  M_save <- lapply(seq_along(S), function(i){
    diff <- (S[[i]] - diag(p))
    priors[[i]] * mat_power(diff, 2)
  }) |>
    Reduce(`+`, x = _)

  M <- M_save - mat_power(M_sir, 2)

  if(regularize){
    if(is.null(lambda)) lambda <- 10e-6
    if(length(lambda) == 1L){
      N <- N + lambda*diag(ncol(N))
    }
  }

  eigs <- geigen::geigen(M, N)
  beta <- as.matrix(with(eigs, vectors[,ncol(vectors):1])[,dims])
  eigvalues <- with(eigs, values[length(values):1])[dims]
  ## Return
  out <- list("ProjectedData" = x %*% beta,
              "ProjectionMatrix" = beta,
              "eigvalues" = eigvalues,
              "M" = M, "N" = N,
              "x" = x, "y" = y,
              "ytype" = ytype, "slices" = slices
  )
  class(out) <- c("sdr")
  out
}


sdr.fit.rsir <- function(object, x, y, ytype, dims, priors = NULL, lambda = 10e-6, ...){
  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }


  data_list <- data_list_fn(x, slices)
  if(is.null(priors)) priors <- priors_fn(data_list = data_list)

  xbarbar <- as.matrix(colMeans(x))
  xbar <- lapply(data_list, function(x) as.matrix(colMeans(x)))
  M <- lapply(seq_along(xbar), function(i){
    diff <- xbar[[i]] - xbarbar
    priors[[i]] * diff %*% t(diff)
  }) |>
    Reduce(`+`, x = _)

  S <- cov(x)
  p <- ncol(S)
  if(length(lambda) == 1){
    N <- S + lambda*diag(p)
    eigs <- geigen::geigen(M, N)
  }

  beta <- as.matrix(with(eigs, vectors[,ncol(vectors):1])[,dims])
  eigvalues <- with(eigs, values[length(values):1])[dims]
  ## Return
  list("ProjectedData" = x %*% beta, "ProjectionMatrix" = beta, "eigvalues" = eigvalues, "M" = M, "N" = N, "slices" = slices)
}

sdr.fit.rsave <- function(object, x, y, ytype, dims, priors = NULL, lambda = 10e-6, ...){
  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }
  p <- ncol(x)
  S_x <- cov(x)
  S_inv_sqrt <- mat_power(S_x, -.5)
  S_sqrt <- mat_power(S_x, .5)
  z <- scale(x, center = TRUE, scale = FALSE) %*% S_inv_sqrt

  data_list <- data_list_fn(z, slices)
  if(is.null(priors)) priors <- priors_fn(data_list = data_list)

  S_zy <- lapply(data_list, cov)
  E <- lapply(seq_along(S_zy), function(i){
    priors[[i]]* mat_power(diag(p) - S_zy[[i]], 2)
  }) |>
    Reduce(`+`, x = _)
  M <- S_sqrt %*% E %*% S_sqrt

  if(length(lambda) == 1){
    N <- S_x + lambda*diag(p)
    eigs <- geigen::geigen(M, N)
  }

  # beta <- S_inv_sqrt %*% eigen(S_inv_sqrt %*% M %*% S_inv_sqrt)$vectors
  eigs <- geigen::geigen(M, N)
  beta <- as.matrix(with(eigs, vectors[,ncol(vectors):1])[,dims])
  eigvalues <- with(eigs, values[length(values):1])[dims]
  ## Return
  list("ProjectedData" = x %*% beta, "ProjectionMatrix" = beta, "eigvalues" = eigvalues, "M" = M, "N" = N, "slices" = slices)
}

sdr.fit.dr <- function(object, x, y, ytype, dims, priors = NULL, ...){
  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }

  p <- ncol(x)
  N <- S_x <- cov(x)
  S_inv_sqrt <- mat_power(S_x, -.5)
  S_sqrt <- mat_power(S_x, .5)
  z <- scale(x, center = TRUE, scale = FALSE) %*% S_inv_sqrt

  data_list <- data_list_fn(z, slices)
  priors <- priors_fn(data_list = data_list)

  S_zy <- lapply(data_list, cov)
  E_zy <- do.call(rbind, lapply(data_list, colMeans))

  M1 <- lapply(seq_along(data_list), function(i) {
    priors[[i]] * mat_power(S_zy[[i]] + E_zy[i, ] %*% t(E_zy[i, ]), 2)
  }) |>
    Reduce(`+`, x = _)

  M2 <- lapply(seq_along(data_list), function(i){
    priors[[i]] * E_zy[i,] %*% t(E_zy[i,])
  }) |>
    Reduce(`+`, x = _)

  E <- 2*M1 + 2*mat_power(M2, 2) + 2*sum(diag(M2))*M2 - 2*diag(p)
  M <- S_sqrt %*% E %*% S_sqrt

  # beta <- S_inv_sqrt %*% eigen(S_inv_sqrt %*% M %*% S_inv_sqrt)$vectors
  eigs <- geigen::geigen(M, N)
  beta <- as.matrix(with(eigs, vectors[,ncol(vectors):1])[,dims])
  eigvalues <- with(eigs, values[length(values):1])[dims]
  ## Return
  out <- list("ProjectedData" = x %*% beta,
              "ProjectionMatrix" = beta,
              "eigvalues" = eigvalues,
              "M" = M, "N" = N,
              "x" = x, "y" = y,
              "ytype" = ytype, "slices" = slices
  )
  class(out) <- c("sdr")
  out
}

sdr.fit.sdrs <- function(object, x, y, ytype, dims, priors = NULL, prec.est = "glasso", ...){
  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }

  classes <- unique(slices)
  data <- lapply(1:length(classes), function(i){
    x[slices == classes[[i]], ]
  })

  n <- sapply(data, nrow)
  n_order <- order(n, decreasing = TRUE)
  n <- n[n_order]
  prior <- lapply(seq_along(n), function(i) n[[i]]/sum(n))
  data <- data[n_order]
  S <- lapply(data, stats::cov)
  xbar <- lapply(data, colMeans)

  # if(is.list(prec.est)) S_inv <- prec.est

  if(is.null(prec.est)) {
    S_inv <- lapply(S, function(x) solve(x, tol = NULL))
  } else if(prec.est == "Haff"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      haff_est(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "Wang"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      wang_est(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "Bodnar"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      p <- ncol(S[[i]])
      HDShOP::InvCovShrinkBGP16(n[[i]], p, diag(p), solve(S[[i]], tol = NULL))$S
    })
  } else if(prec.est == "MRY"){
    S_inv <- SCPME_qda(x, slices, ...)
    lam <- sapply(seq_along(S_inv), function(i) S_inv[[i]]$Tuning[2])
    S_inv <- lapply(seq_along(S_inv), function(i) S_inv[[i]]$Omega)
  } else if(prec.est == "glasso"){
    glasso_wrapper <- function(data, S, lam = NULL, ...){
      out <- vector(mode = "list", length = length(S))

      S_biased <- lapply(seq_along(S), function(i) S[[i]]*(n[[i]] - 1)/n[[i]])
      # if(is.null(lam)){
      #   lam <- lapply(S_biased, function(x){
      #     eig_vals <- eigen(x)$values
      #     max_off_diag <- max(abs(x[upper.tri(x)]))
      #     exp(seq(
      #       log(min(eig_vals)),
      #       log(min(max(eig_vals), max_off_diag)),
      #       length.out = 10)
      #     )
      #   })
      # }

      for (i in seq_along(S)) {

        out[[i]] <- CVglasso::CVglasso(data[[i]], S = S_biased[[i]],
                                       trace = "none", crit.cv = "loglik", lam = lam[[i]])
        # if(is.null(lam)){
        #   min.ratio <- 0.01
        #   while(temp$Lambdas[1] == temp$Tuning[[2]] && min.ratio >= 0.000625){
        #     min.ratio <- min.ratio/2
        #     temp <- CVglasso::CVglasso(x, trace = "none", crit.cv = "loglik", lam.min.ratio = min.ratio)
        #   }
        # }
        # out[[i]] <- temp
      }
      out
    }
    glasso_out <- glasso_wrapper(data, S, ...)
    # S <- lapply(glasso_out, function(x) x$Sigma)
    S_inv <- lapply(glasso_out, function(x) x$Omega)
    lam <- sapply(glasso_out, function(x) x$Tuning[[2]])
  } else if(prec.est == "pseudo"){
    S_inv <- lapply(S, MASS::ginv)
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
  M_svd <- svd(M)
  beta <- as.matrix(M_svd$u[,dims])
  eigvalues <- (M_svd$d^2)[dims]
  #Return Projection Matrix

  ## Return
  out <- list("ProjectedData" = x %*% beta,
              "ProjectionMatrix" = beta,
              "eigvalues" = eigvalues,
              "M" = M, "N" = diag(ncol(M)),
              "x" = x, "y" = y,
              "ytype" = ytype, "slices" = slices
  )
  class(out) <- c("sdr")

  if(!is.null(prec.est)){
    if(prec.est == "SCPME" | prec.est == "glasso") out$lam <-  lam
  }
  out
}

sdr.fit.sdrs2 <- function(object, x, y, ytype, dims, priors = NULL, prec.est = "glasso", ...){
  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }

  classes <- unique(slices)
  data <- lapply(1:length(classes), function(i){
    x[slices == classes[[i]], ]
  })

  n <- sapply(data, nrow)
  n_order <- order(n, decreasing = TRUE)
  n <- n[n_order]
  prior <- lapply(seq_along(n), function(i) n[[i]]/sum(n))
  data <- data[n_order]
  S <- lapply(data, stats::cov)
  # S_common <- cov(x)
  xbar <- lapply(data, colMeans)
  xbarbar <- colMeans(x)

  # if(is.list(prec.est)) S_inv <- prec.est

  if(is.null(prec.est)) {
    S_inv <- lapply(S, function(x) solve(x, tol = NULL))
  } else if(prec.est == "Haff"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      Haffshrinkage(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "Wang"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      wang_est(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "Bodnar"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      p <- ncol(S[[i]])
      HDShOP::InvCovShrinkBGP16(n[[i]], p, diag(p), solve(S[[i]], tol = NULL))$S
    })
  } else if(prec.est == "MRY"){
    S_inv <- SCPME_qda(x, slices, ...)
    lam <- sapply(seq_along(S_inv), function(i) S_inv[[i]]$Tuning[2])
    S_inv <- lapply(seq_along(S_inv), function(i) S_inv[[i]]$Omega)
  } else if(prec.est == "glasso"){

    glasso_wrapper <- function(data, S, lam = NULL, ...){
      out <- vector(mode = "list", length = length(S))

      S_biased <- lapply(seq_along(S), function(i) S[[i]]*(n[[i]] - 1)/n[[i]])
      # if(is.null(lam)){
      #   lam <- lapply(S_biased, function(x){
      #     eig_vals <- eigen(x)$values
      #     eig_zero <- sapply(eig_vals, is.zero)
      #     eig_vals <- eig_vals[!eig_zero]
      #     max_off_diag <- max(abs(x[upper.tri(x)]))
      #     exp(seq(
      #       log(min(eig_vals)),
      #       log(min(max(eig_vals), max_off_diag)),
      #       length.out = 10)
      #     )
      #   })
      # }

      for (i in seq_along(S)) {

        out[[i]] <- CVglasso::CVglasso(data[[i]], S = S_biased[[i]],
                                       trace = "none", crit.cv = "loglik", lam = lam[[i]])
        # if(is.null(lam)){
        #   min.ratio <- 0.01
        #   while(temp$Lambdas[1] == temp$Tuning[[2]] && min.ratio >= 0.000625){
        #     min.ratio <- min.ratio/2
        #     temp <- CVglasso::CVglasso(x, trace = "none", crit.cv = "loglik", lam.min.ratio = min.ratio)
        #   }
        # }
        # out[[i]] <- temp
      }
      out
    }
    glasso_out <- glasso_wrapper(data, S, ...)
    # S <- lapply(glasso_out, function(x) x$Sigma)
    S_inv <- lapply(glasso_out, function(x) x$Omega)
    lam <- sapply(glasso_out, function(x) x$Tuning[[2]])

  }

  # S <- lapply(S_inv, function(x) solve(x))

  projectedMeanDiffs <- do.call(cbind,
                                lapply(1:length(xbar), function(i) {
                                  S_inv[[i]] %*% (xbar[[i]] - xbarbar)
                                }))

  Sdiffs <- do.call(cbind,
                    lapply(2:length(S), function(i) {
                      S[[i]] - S[[1]]
                    }))

  Sdiff_svd <- svd(Sdiffs)

  top <- lik(Sdiff_svd$d^2)$max

  # out <- core(Sigmas = S, ns = n, numdir.test = TRUE, numdir = ncol(S[[1]]))
  # d <- which.min(out$bic)
  # Sdiffs <- out$Gammahat[[d]]

  # Sdiifs <- with(Sdiff_svd,
  #                as.matrix(u[,1:top]) %*% as.matrix(diag(d))[1:top, 1:top])

  Sdiffs <- Sdiff_svd$u[,1:top]

  M <- cbind(apply(projectedMeanDiffs, 2, normalize), Sdiffs)
  M_svd <- svd(M, nu = nrow(M))

  beta <- as.matrix(M_svd$u[,dims])
  eigvalues <- (M_svd$d^2)[dims]

  ## Return
  out <- list("ProjectedData" = x %*% beta, "ProjectionMatrix" = beta, "eigvalues" = eigvalues, "M" = M, "N" = diag(nrow(M)), "slices" = slices)
  # out <- list("ProjectionMatrix" = M_svd$u, "eigvalues" = M_svd$d^2, "M" = M, "N" = diag(nrow(M)), "slices" = slices)

  if(!is.null(prec.est)){
    if(prec.est == "SCPME" | prec.est == "glasso") out$lam <-  lam
  }
  out
}


sdr.fit.sdrs3 <- function(object, x, y, ytype, dims, priors = NULL, ...){

  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }

  classes <- unique(slices)
  data <- lapply(1:length(classes), function(i){
    x[slices == classes[[i]], ]
  })

  n <- sapply(data, nrow)
  n_order <- order(n, decreasing = TRUE)
  n <- n[n_order]
  prior <- lapply(seq_along(n), function(i) n[[i]]/sum(n))
  data <- data[n_order]
  S_x <- cov(x)
  S_inv_sqrt <- mat_power(S_x, power = -1/2)
  S <- lapply(data, stats::cov)

  # S_common <- cov(x)
  xbar <- lapply(data, function(x)as.matrix(colMeans(x)))
  xbarbar <- as.matrix(colMeans(x))

  # mean_diffs <- Reduce(`+`,
  #                      lapply(1:length(xbar), function(i) {
  #                        diff <- xbar[[i]] - xbarbar
  #                        prior[[i]] * (diff %*% t(diff))
  #                        })
  #                      )
  mean_diffs <- do.call(cbind,
                       lapply(1:length(xbar), function(i) {
                         S_inv_sqrt %*% (xbar[[i]] - xbarbar)
                         })
                       )


  Sdiffs <- do.call(cbind,
                    lapply(2:length(S), function(i) {
                      S[[i]] - S[[1]]
                    }))

  Sdiff_svd <- corpcor::fast.svd(Sdiffs)

  top <- lik(Sdiff_svd$d^2)$max

  # out <- core(Sigmas = S, ns = n, numdir.test = TRUE, numdir = ncol(S[[1]]))
  # d <- which.min(out$bic)
  # Sdiffs <- out$Gammahat[[d]]

  # Sdiifs <- with(Sdiff_svd,
  #                as.matrix(u[,1:top]) %*% as.matrix(diag(d))[1:top, 1:top])

  Sdiffs <- Sdiff_svd$u[,1:top]

  M <- cbind(mean_diffs, Sdiffs)
  M_svd <- svd(M, nu = nrow(M))

  beta <- as.matrix(M_svd$u[,dims])
  eigvalues <- (M_svd$d^2)[dims]
  #Return Projection Matrix

  out <- list("ProjectedData" = x %*% beta,
              "ProjectionMatrix" = beta,
              "eigvalues" = eigvalues,
              "M" = M, "N" = diag(ncol(M)),
              "x" = x, "y" = y,
              "ytype" = ytype, "slices" = slices
  )
  class(out) <- c("sdr")
  out
}

sdr.fit.pfc <- function(object, x, y, ytype, dims, ...){
  if(ytype == "categorical"){
    slices <- y
  } else {
    slices <- make_slices(y, ...)
  }

  fy <- model.matrix(~ slices - 1)[, 1:(length(unique(slices)) - 1)]
  fy <- as.matrix(scale(fy, center = TRUE, scale = FALSE))

  xc <- scale(x, TRUE, FALSE)
  P_F <- as.matrix(fy %*% solve(t(fy) %*% fy) %*% t(fy))
  N <- S <- cov(xc)
  M <- S_fit <- cov(P_F %*% x)

  eigs <- geigen::geigen(M, N)
  beta <- as.matrix(with(eigs, vectors[,ncol(vectors):1])[,dims])
  eigvalues <- with(eigs, values[length(values):1])[dims]
  out <- list("ProjectedData" = x %*% beta,
              "ProjectionMatrix" = beta,
              "eigvalues" = eigvalues,
              "M" = M, "N" = diag(ncol(M)),
              "x" = x, "y" = y,
              "ytype" = ytype, "slices" = slices
  )
  class(out) <- c("sdr")
  out
}

sdr.fit.default <- function(method, ...) {
  stop(paste0("Unrecognized 'method = ", method, "'. See ?sdr for available methods."))
}
