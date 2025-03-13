core <- function(X, y, Sigmas = NULL, ns = NULL, numdir = 2, numdir.test = FALSE, ...) {
  mf <- match.call()

  if (!is.null(Sigmas)) {
    if (is.null(ns))
      stop("Number of observations per class must be provided")
    p <- dim(Sigmas[[1]])[1]
    hlevels <- length(Sigmas)
  }

  if (is.null(Sigmas) & is.null(ns)) {
    if (is.factor(y))
      y <- as.integer(y)
    hlevels <- length(unique(y))
    p <- ncol(X)
    Sigmas <- list()
    ns <- vector(length = hlevels)
    for (h in 1:hlevels) {
      ns[h] <- sum(y == h)
      Sigmas[[h]] <- cov(as.matrix(X[y == h, ], ncol = p))
    }
  }

  n <- sum(ns)
  ff <- ns / n
  d <- numdir

  orthonorm <- function(u) {
    if (is.null(u))
      return(NULL)
    if (!(is.matrix(u)))
      u <- as.matrix(u)
    dd <- dim(u)
    n <- dd[1]
    p <- dd[2]
    if (prod(abs(La.svd(u)$d) > 1e-08) == 0)
      stop("collinears vectors in orthonorm")
    if (n < p) {
      warning("There are too much vectors to orthonormalize in orthonorm.")
      u <- as.matrix(u[, 1:p])
      n <- p
    }
    v <- u
    if (p > 1) {
      for (i in 2:p) {
        coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)])) / diag(crossprod(v[, 1:(i - 1)]))
        v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = n) %*% matrix(coef.proj, nrow = i - 1)
      }
    }
    coef.proj <- 1 / sqrt(diag(crossprod(v)))
    return(t(t(v) * coef.proj))
  }

  projection <- function(alpha, Sigma)
  {
    alpha %*% solve(t(alpha) %*% Sigma %*% alpha) %*%
      t(alpha) %*% Sigma
  }

  one.core <- function(d) {
    if (d == 0) {
      Gamma.hat <- NULL
      Sigma.hat <- matrix(0, p, p)
      for (g in 1:hlevels) {
        Sigma.hat <- Sigma.hat + ff[g] * Sigmas[[g]]
      }
      term0 <- 0
      term1 <- n / 2 * log(det(Sigma.hat))
      loglik <- term0 - term1
    } else if (d == p) {
      Gamma.hat <- diag(p)
      Sigma.hat <- matrix(0, p, p)
      for (g in 1:hlevels) {
        Sigma.hat <- Sigma.hat + ff[g] * Sigmas[[g]]
      }
      term0 <- 0
      term1 <- n / 2 * log(det(Sigma.hat))
      term2 <- n / 2 * log(det(Sigma.hat))
      term3 <- 0
      for (g in 1:hlevels) {
        term3 <- term3 + ns[g] / 2 * log(det(Sigmas[[g]]))
      }
      loglik <- Re(term0 - term1 + term2 - term3)
    } else {
      objfun <- function(W) {
        Q <- W$Qt
        d <- W$dim[1]
        p <- W$dim[2]
        Sigmas <- W$Sigmas
        n <- sum(W$ns)
        U <- matrix(Q[, 1:d], ncol = d)
        V <- matrix(Q[, (d + 1):p], ncol = (p - d))
        Sigma.hat <- matrix(0, p, p)
        for (g in 1:hlevels) {
          Sigma.hat <- Sigma.hat + ff[g] * Sigmas[[g]]
        }
        Ps <- projection(U, diag(p))
        term0 <- 0
        term1 <- n / 2 * log(det(Sigma.hat))
        term2 <- n / 2 * log(det(t(U) %*% Sigma.hat %*% U))
        term3 <- 0
        for (g in 1:hlevels) {
          term3 <- term3 + ns[g] / 2 * log(det(t(U) %*% Sigmas[[g]] %*% U))
        }
        value <- Re(term0 - term1 + term2 - term3)
        return(list(value = value))
      }

      W <- list(dim = c(d, p), Sigmas = Sigmas, ns = ns)
      grassmann <- GrassmannOptim(objfun, W, ...)  # pass objfun here
      Gamma.hat <- matrix(grassmann$Qt[, 1:d], ncol = d)
      loglik <- tail(grassmann$fvalues, n = 1)
    }

    Sigma.hat <- matrix(0, p, p)
    for (g in 1:hlevels) {
      Sigma.hat <- Sigma.hat + ff[g] * Sigmas[[g]]
    }
    if (d != 0) {
      Ps.hat <- projection(Gamma.hat, Sigma.hat)
    }
    Sigmas.hat <- list()
    for (g in 1:hlevels) {
      if (d == 0) {
        Sigmas.hat[[g]] <- Sigma.hat
      } else {
        Sigmas.hat[[g]] <- Sigma.hat + t(Ps.hat) %*% (Sigmas[[g]] - Sigma.hat) %*% Ps.hat
      }
    }

    numpar <- p * (p + 1) / 2 + d * (p - d) + (hlevels - 1) * d * (d + 1) / 2
    aic <- -2 * loglik + 2 * numpar
    bic <- -2 * loglik + log(n) * numpar
    return(list(Gammahat = Gamma.hat, Sigmahat = Sigma.hat, Sigmashat = Sigmas.hat, loglik = loglik, numpar = numpar, aic = aic, bic = bic))
  }

  if (!numdir.test) {
    fit <- one.core(d)
    ans <- list(Gammahat = fit$Gammahat, Sigmahat = fit$Sigmahat, Sigmashat = fit$Sigmashat, loglik = fit$loglik, aic = fit$aic, bic = fit$bic, numpar = fit$numpar, numdir = d, model = "core", call = match.call(expand.dots = TRUE), numdir.test = numdir.test)
    class(ans) <- "core"
    return(invisible(ans))
  }

  aic <- bic <- numpar <- loglik <- vector(length = d)
  Gammahat <- Sigmahat <- Sigmashat <- list()
  for (i in 1:d) {
    if (!is.null(mf$verbose))
      cat("Running CORE for numdir =", i, "\n")
    fit <- one.core(i)
    Gammahat[[i]] <- fit$Gammahat
    Sigmahat[[i]] <- fit$Sigmahat
    Sigmashat[[i]] <- fit$Sigmashat
    loglik[i] <- fit$loglik
    numpar[i] <- fit$numpar
    aic[i] <- fit$aic
    bic[i] <- fit$bic
  }

  ans <- list(Gammahat = Gammahat, Sigmahat = Sigmahat, Sigmashat = Sigmashat, loglik = loglik, aic = aic, bic = bic, numpar = numpar, numdir = d, model = "core", call = match.call(expand.dots = TRUE), numdir.test = numdir.test)
  class(ans) <- "core"
  return(invisible(ans))
}
