qdap <- function(x, y, prior){
  vecs <- list(x, y)
  s2 <- sapply(vecs, var)
  ords <- sort(s2, decreasing = TRUE, index.return = TRUE)$ix
  # ords <- order(unlist(s2), decreasing = TRUE)
  s2 <- s2[ords]
  # s2 <- as.list(unlist(s2)[ords])
  # s <- lapply(s2, sqrt)
  s <- sqrt(s2)
  prior <- prior[ords]
  mdiff <- mean(vecs[[ords[1]]]) - mean(vecs[[ords[2]]])
  s2diff <- Reduce(`-`, s2)
  sqrt_delt <- sqrt((-1*mdiff)^2 + s2diff*log(s2[[1]]/s2[[2]]))
  sdelt <- lapply(s, function(x) x*sqrt_delt)
  smdiff <- lapply(s, function(x) x*mdiff)
  # Return
  c("qdap" =
      prior[[1]] +
      prior[[1]]*pnorm((smdiff[[2]] - sdelt[[1]])/s2diff) -
      prior[[1]]*pnorm((smdiff[[2]] + sdelt[[1]])/s2diff) +
      prior[[2]]*pnorm((smdiff[[1]] + sdelt[[2]])/s2diff) -
      prior[[2]]*pnorm((smdiff[[1]] - sdelt[[2]])/s2diff)
    )

}

# dim_order_cv <- function(x, grouping, cv_method = "t", ndims = NULL, folds = 10, nsims = 1, npca = NULL, lambdas = NULL, ...){
#   grouping <- droplevels(as.factor(grouping))
#   classes <- unique(grouping)
#   if(length(classes) > 2) stop("dim_order only currently supports binary classes")
#   if(is.null(ndims)) ndims <- ncol(x)
#   df <- bind_cols("class" = grouping, x)
#
#
#   parallel::mclapply(1:nsims, function(sim){
#     set.seed(sim)
#     splits <- rsample::vfold_cv(df, strata = class, v = folds)$splits
#     train <- lapply(splits, rsample::training)
#     test <- lapply(splits, rsample::testing)
#
#     temp <- matrix(NA, nrow = length(train), ncol = npca)
#
#     for(i in seq_along(train)){
#       tr <- train[[i]]
#       ts <- test[[i]]
#
#       ## PCA
#       if(!is.null(npca)){
#         ## PCA decomp
#         pca_result <- prcomp(tr[,-1], center = TRUE, scale. = TRUE)
#         ## PCA Porjection Matrix
#         pca_proj <- pca_result$rotation[,1:npca]
#         ## Test scaled and projected
#         train_mean <- colMeans(tr[,-1])
#         train_sd <- apply(tr[,-1], 2, sd)
#         ts <- bind_cols("class" = ts$class,
#                         as.matrix(scale(ts[,-1], center = train_mean, scale = train_sd)) %*% pca_proj)
#         ## Train projected
#         tr <- bind_cols("class" = tr$class, as.matrix(tr[,-1]) %*% pca_proj)
#         ndims <- npca
#       }
#
#       mod <- hldr::hldr(class ~., data = tr, method = "SYS", prec.est = "SCPME", lam = lambdas)
#       temp[i,1:ndims] <- dim_order(x = ts[,-1],
#                             grouping = ts$class, method = cv_method,
#                 ProjectionMatrix = mod$ProjectionMatrix)$dim_crit
#
#     }
#     temp
#   }, mc.cores = parallel::detectCores() - 1) |>
#     do.call(rbind, args = _) |>
#     colMeans()
# }

dim_order <- function(x, grouping, ProjectionMatrix, method = "t",
                      ndims = NULL, ProjectedData = NULL, prior = NULL, nboot = NULL, object = NULL, ...){

  method <- tolower(method)
  if(!(method %in% c("t", "qdap", "f"))) stop(paste0("Unrecognized method = ", method, '. dim_order only supports methods "T", "QDAP", or "F".'))
  if(!is.null(object)){
    if(class(object) != "sdr") stop('dim_order only supports object class of "sdr"')
    ProjectionMatrix <- object$ProjectionMatrix
    ProjectedData <- object$ProjectedData
    grouping <- object$slices
  }

  if(!is.factor(grouping))  grouping <- droplevels(as.factor(grouping))
  classes <- unique(grouping)
  if((method == "t" || method == "qdap") && length(classes) > 2){
    stop(paste0("method = ", method, ' only supports binary classes. Maybe try method = "F"?'))
  }

  if(is.null(object)){
    if(missing(grouping)) stop("If object = NULL, must supply grouping argument")
    if(is.null(ProjectedData)){
      if(missing(ProjectionMatrix) || missing(x)){
        stop("If ProjectedData = NULL, must supply x, grouping, and ProjectionMatrix arguments")
      }
      if(!is.matrix(x)) x <- as.matrix(x)
      if(!is.matrix(ProjectionMatrix)) ProjectedMatrix <- as.matrix(ProjectionMatrix)
      ProjectedData <- x %*% ProjectionMatrix
    }
  }
  if(is.null(ndims)) ndims <- ncol(ProjectedData)

  if(!is.null(nboot)){

    # df <- cbind("class" = grouping, ProjectedData) |> as.data.frame()
    # splits <- (df |> rsample::vfold_cv(strata = class, v = folds))$splits
    # train <- lapply(splits, rsample::training)
    # test <- lapply(splits, rsample::testing)

    crit <- lapply(seq_len(nboot), function(i){
      # train <- train[[i]]
      # proj_list <- train |>
      #   group_by(class) |>
      #   group_split() |>
      #   lapply(function(x) as.matrix(x[,-1]))
      #
      # proj_list <- data_list_fn(train[[i]][,-1], train)

      idx <- sample(1:nrow(ProjectedData), nrow(ProjectedData), replace = TRUE)
      dim_order.fit(ProjectedData[idx, ], grouping[idx], method, ...)
      # do.call(
      #   cbind,
      #   lapply(1:ncol(train[,-1]), function(r){
      #     abs(t.test(proj_list[[1]][,r], proj_list[[2]][,r],
      #                alternative = "two.sided",
      #                var.equal = FALSE)$statistic)
      #   })
      # )
    }) |>
      do.call(rbind, args = _) |>
      colMeans()

  } else {
  crit <- dim_order.fit(ProjectedData, grouping, method, ...)
  }
  if(method == "qdap"){
    dcr <- FALSE
  } else {
    dcr <- TRUE
  }
  topdims <- sort(crit, decreasing = dcr, index.return = TRUE)

  list("dims" = topdims$ix[1:ndims],
       "dim_criteria" = crit[1:ndims],
       "ProjectionMatrix" = ProjectionMatrix[,topdims$ix],
       "ProjectedData" = ProjectedData[,topdims$ix])
}

dim_order.fit <- function(ProjectedData, grouping, method, ...){
  if(method == "t"){
    proj_list <- data_list_fn(ProjectedData, grouping)
    crit <-
      sapply(1:ncol(ProjectedData), function(r){
        abs(t.test(proj_list[[1]][,r], proj_list[[2]][,r],
                   alternative = "two.sided",
                   var.equal = FALSE)$statistic)
      })
    # topdims <- sort(crit, decreasing = TRUE, index.return = TRUE)
  }
  if(method == "qdap"){
    proj_list <- data_list_fn(ProjectedData, grouping)
    if(is.null(prior)) prior <- as.vector(table(grouping))/length(grouping)
    crit <- sapply(1:ncol(ProjectedData), function(r){
      qdap(proj_list[[1]][,r], proj_list[[2]][,r], prior)
    })
    # topdims <- sort(crit, decreasing = FALSE, index.return = TRUE)
  }
  if(method == "f"){
    crit <-
      sapply(1:ncol(ProjectedData), function(r){
        df <- data.frame("proj_data_r" = ProjectedData[,r],
                         "grouping" = grouping)
        oneway.test(proj_data_r ~ grouping,
                    data = df, var.equal = FALSE)$statistic
      })
    # topdims <- sort(crit, decreasing = TRUE, index.return = TRUE)
  }
  crit
}
