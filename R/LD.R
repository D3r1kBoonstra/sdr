# # Helper Functions --------------------------------------------------------
# logd <- function(matrix){
#   matSVD <- svd(matrix)
#   u <- matSVD$u
#   v <- matSVD$v
#   u %*% diag(log(matSVD$d), nrow = length(matSVD$d)) %*% t(v)
# }
#
# matInvSqrt <- function (matrix) {
#   svdMat <- svd(matrix)
#   U <- svdMat$u
#   rootDInv <- sqrt(solve(diag(svdMat$d, nrow = length(svdMat$d)), tol = NULL))
#   U %*% rootDInv %*% t(U)
# }
#
#
#
#
# # LD Function -------------------------------------------------------------
# LD <- function(x, grouping, dims, prior = NULL, ...){
#
#   grouping <- droplevels(as.factor(grouping))
#   classes <- unique(grouping)
#   data <- lapply(1:length(classes), function(i){
#     x[grouping == classes[[i]], ]
#   })
#
#   prior <- prior(data)
#   if(is.null(prior)) prior <- as.list((table(grouping))/length(grouping))
#   xbar <- lapply(data, colMeans)
#   B <- S_B(prior, xbar)
#   S <- S_W(prior, data)
#   covs <- lapply(data, cov)
#
#   combns <- combn(length(data), 2, simplify = FALSE)
#   mld_diff <- lapply(combns, function(x){
#     (xbar[[x[1]]] - xbar[[x[2]]]) %*% t(xbar[[x[1]]] - xbar[[x[2]]])
#   })
#
#   mld_pie <- lapply(combns, function(x){
#     c(prior[[x[1]]] / (prior[[x[1]]] + prior[[x[2]]]),
#       prior[[x[2]]] / (prior[[x[1]]] + prior[[x[2]]]))
#   })
#
#   Sij <- lapply(1:length(combns), function(x){
#     combns <- combns[[x]]
#     mld_pie[[x]][1] * covs[[combns[1]]] +
#       mld_pie[[x]][2] * covs[[combns[2]]]
#   })
#
#   mld_fun  <- lapply(1:length(combns), function(x){
#     combns <- combns[[x]]
#     someRootInv <- matInvSqrt(
#       matInvSqrt(S) %*%
#         Sij[[x]] %*%
#         matInvSqrt(S)
#     )
#
#     prior[[combns[1]]] * prior[[combns[2]]] * solve(S, tol = NULL) %*% solve(matInvSqrt(S), tol = NULL) %*%
#       (
#         someRootInv %*%
#           matInvSqrt(S) %*%
#           mld_diff[[x]] %*%
#           matInvSqrt(S) %*%
#           someRootInv +
#           (1 / (mld_pie[[x]][1] * mld_pie[[x]][2])) * logd(matInvSqrt(S) %*%
#                                                              Sij[[x]] %*%
#                                                              matInvSqrt(S)) -
#           mld_pie[[x]][1] *  logd(matInvSqrt(S) %*%
#                                     covs[[combns[1]]] %*%
#                                     matInvSqrt(S)) -
#           mld_pie[[x]][2] *  logd(matInvSqrt(S) %*%
#                                     covs[[combns[2]]] %*%
#                                     matInvSqrt(S))
#
#       ) %*%
#       solve(matInvSqrt(S), tol = NULL)
#   })
#
#   M <- Reduce(`+`, mld_fun)
#   U <- as.matrix(svd(M)$u)
#   #Return Projection Matrix
#   list("M" = M, "ProjectionMatrix" = t(U[,dims]))
# }
#
#
# SAVE <- function(x, grouping, dims, prior = NULL, ...){
#   grouping <- droplevels(as.factor(grouping))
#   classes <- unique(grouping)
#   data <- lapply(1:length(classes), function(i){
#     x[grouping == classes[[i]], ]
#   })
#
#   prior <- prior(data)
#   if(is.null(prior)) prior <- as.list((table(grouping))/length(grouping))
#   xbar <- lapply(data, colMeans)
#   B <- S_B(prior, xbar)
#   S <- S_W(prior, data)
#   covs <- lapply(data, cov)
#
#   S_hatGamma <- (1 / length(covs)) *
#     Reduce(`+`,
#            lapply(covs, function(x){
#              (x - S) %*%
#                solve(B + S, tol = NULL) %*%
#                (x - S)
#            })
#     )
#
#   rootGammaInv <- matInvSqrt(B + S)
#
#   M <- (rootGammaInv %*% B %*% rootGammaInv) %*%
#     (rootGammaInv %*% B %*% rootGammaInv) +
#     (rootGammaInv %*% S_hatGamma %*% rootGammaInv)
#   U <- rootGammaInv %*% as.matrix(svd(M)$u)
#   list("M" = M, "ProjectionMatrix" = t(U[,dims]))
# }
#
# SIR <- function(x, grouping, dims, prior = NULL, ...){
#   grouping <- droplevels(as.factor(grouping))
#   classes <- unique(grouping)
#   data <- lapply(1:length(classes), function(i){
#     x[grouping == classes[[i]], ]
#   })
#
#   prior <- prior(data)
#   if(is.null(prior)) prior <- as.list((table(grouping))/length(grouping))
#   xbar <- lapply(data, colMeans)
#   B <- S_B(prior, xbar)
#   S <- S_W(prior, data)
#   rootGammaInv <- matInvSqrt(B + S)
#   M <- rootGammaInv %*% B %*% rootGammaInv
#   U <- rootGammaInv %*% do.call(svd, list(M))$u
#   list("M" = M, "ProjectionMatrix" = t(U[,dims]))
#
# }
