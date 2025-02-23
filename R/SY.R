SY <- function(x, grouping, dims, ...){
  classes <- unique(grouping)
  data <- lapply(1:length(classes), function(i){
    x[grouping == classes[[i]], ]
  })

  n <- lapply(data, nrow)
  S <- lapply(data, stats::cov)
  xbar <- lapply(data, colMeans)

  S_inv <- lapply(seq_along(S), function(i, ...){
    solve(S[[i]], tol = NULL)
  })

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
  list("M" = M, "ProjectionMatrix" = t(U[,dims]))
}
