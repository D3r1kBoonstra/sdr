sdr <- function(x, y, method = "sdrs", ytype = "categorical", dims = NULL, dim.order = NULL, cv.folds = NULL, ...){

  if(!is.matrix(x)) x <- as.matrix(x)
  if(tolower(ytype) == "categorical") slices <- y
  if(!is.factor(slices)) slices <- as.factor(slices)
  slices <- droplevels(slices)
  if(is.null(dims)) dims <- 1:ncol(x)


  out <- sdr.fit(x = x, slices = slices, method = method, dims = dims, dim.order = dim.order, ...)

  if(!is.null(cv.folds)){
    df <- cbind("class" = slices, as.data.frame(x))
    splits <- rsample::vfold_cv(df, strata = class, v = cv.folds)$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    cv.errs <- matrix(nrow = cv.folds, ncol = length(out$dims))
    for(i in 1:cv.folds){
      fit <- sdr.fit(x = as.matrix(train[[i]][,-1]), slices = train[[i]]$class, method = method, dims = out$dims, dim.order = NULL, lam = out$lam, ...)
      cv.errs[i,] <- sapply(1:length(out$dims), function(j){
        pred <- predict(fit, newdata = test[[i]][,-1], dims = out$dims[1:j], ...)$class
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
    out$ProjectionMatrix <- out$ProjectionMatrix[,1:min_dims]
    out$ProjectedData <- as.matrix(out$ProjectedData[,1:min_dims])
    out$cv.errs <- mean_errs
  }


  #### Model output
  out$call <- match.call()
  out$method <- method
  out$dim_order_method <- dim.order
  # out$dims <- if(!is.null(dimselect)) out$dims
  # out$slices <- slices
  ## Return
  out
}
