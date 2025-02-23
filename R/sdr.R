sdr <- function(x, grouping, method = "SDRS", dims = NULL, dimselect = NULL, cv.folds = NULL, ...){

  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.factor(grouping)) grouping <- as.factor(grouping)
  grouping <- droplevels(grouping)
  if(is.null(dims)) dims <- 1:ncol(x)


  out <- sdr.fit(x = x, grouping = grouping, method = method, dims = dims, dimselect = dimselect, ...)

  if(!is.null(cv.folds)){
    df <- cbind("class" = grouping, as.data.frame(x))
    splits <- rsample::vfold_cv(df, strata = class, v = cv.folds)$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    cv.errs <- matrix(nrow = cv.folds, ncol = length(out$dims))
    for(i in 1:cv.folds){
      fit <- sdr.fit(x = train[[i]][,-1], grouping = train[[i]]$class, method = method, dims = out$dims, dimselect = NULL, lam = out$lam, ...)
      cv.errs[i,] <- sapply(1:length(out$dims), function(j){
        pred <- predict(fit, newdata = test[[i]][,-1], dims = 1:j)$class
        mean(pred != test[[i]]$class)
      })
    }
    mean_errs <- colMeans(cv.errs)
    min_dims <- which(mean_errs == min(mean_errs))
    if(length(min_dims) > 1){
      min_sd <- which.min(apply(cv.errs[,min_dims], 2, sd))
      min_dims <- min_dims[min_sd]
    }
    out$dims <- out$dims[1:min_dims]
    out$ProjectionMatrix <- out$ProjectionMatrix[1:min_dims,]
    out$ProjectedData <- as.matrix(out$ProjectedData[,1:min_dims])
    out$cv.errs <- mean_errs
  }


  #### Model output
  out$call <- match.call()
  out$method <- method
  out$dimselect <- dimselect
  # out$dims <- if(!is.null(dimselect)) out$dims
  out$grouping <- grouping
  ## Return
  out
}
