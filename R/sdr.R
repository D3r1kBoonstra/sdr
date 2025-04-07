sdr <- function(x, y, method = "sdrs", ytype = "guess", nslices = 5, ndims = NULL, dims = NULL, dim.order = NULL, cv.folds = NULL, ...){

  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.numeric(x)) x <- apply(x, 2, as.numeric)

  if (!is.null(colnames(x))) {
    var_names <- colnames(x)
  } else {
    var_names <- NULL
  }
  dimnames(x) <- NULL

  ytype_o <- ytype
  ytype <- tolower(ytype)
  if (!(ytype %in% c("guess", "categorical", "continuous"))) {
    warning(paste0("Invalid ytype = '", ytype_o, "'. Defaulting to ytype = 'guess'"))
    ytype <- "guess"
  }

  # Automatic guessing if ytype is "guess"
  if (ytype == "guess") {
    if (is.factor(y) || is.character(y) || all(table(y)/length(y) > 0.10)) {
      ytype <- "categorical"
    } else {
      ytype <- "continuous"
    }
  }

  # Processing based on ytype
  if (ytype == "categorical") {
    if (!is.factor(y)) y <- as.factor(y)
    y <- droplevels(y)
  } else if (ytype == "continuous") {
    if (!is.numeric(y)) y <- as.numeric(y)
  }


  if(is.null(ndims)) ndims <- ncol(x)
  if(is.null(dims)){
    if(is.null(dim.order)){
      dims <- 1:ndims
    } else {
      dims <- 1:ncol(x)
    }
  }

  extra_args <- list(...)
  out <- sdr.fit(method, x = x, y = y, ytype = ytype, dims = dims, ...)
  out$method <- method
  out$dims <- dims
  out$call <- match.call()
  out$extra_args <- extra_args


  if(!is.null(dim.order)) {
    dim_order_out <- dim_order(object = out, method = dim.order, ndims = ndims, ...)
    out$ProjectionMatrix <- dim_order_out$ProjectionMatrix
    out$dim_order_method <- dim.order
    out$dims <- dim_order_out$dims
    out$eigvalues <- out$eigvalues[dims]
    out$ProjectedData <- dim_order_out$ProjectedData
    out$dim_criteria <- dim_order_out$dim_criteria
    if(!is.null(dim_order_out$boot_crits)) out$boot_crits <- dim_order_out$boot_crits
  }
  # else{
  #   dims <- dims[1:ndims]
  #   out$eigvalues <- out$eigvalues[dims]
  #   out$dims <- dims
  #   out$ProjectionMatrix <- as.matrix(out$ProjectionMatrix[,dims])
  #   out$ProjectedData <- out$ProjectedData[,dims]
  # }

  if(!is.null(cv.folds)){
    df <- cbind("y" = y, as.data.frame(x))
    strata <- NULL
    if(ytype == "categorical") strata <- y
    splits <- rsample::vfold_cv(df, strata = strata, v = cv.folds)$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    cv.errs <- matrix(nrow = cv.folds, ncol = length(out$dims))
    for(i in 1:cv.folds){
      fit <- sdr.fit(method, x = as.matrix(train[[i]][,-1]), y = train[[i]]$y, ytype = ytype, dims = out$dims, ...)
  #     fit <- sdr.fit(x = as.matrix(train[[i]][,-1]), slices = train[[i]]$class, method = method, dims = out$dims, dim.order = NULL, lam = out$lam, ...)
      cv.errs[i,] <- sapply(1:length(out$dims), function(r){
        pred <- predict(fit, newdata = test[[i]][,-1], ndims = r, ...)
        if(ytype == "categorical") {
          mean(pred[[1]] != test[[i]]$y)
        } else {
          Metrics::rmse(test[[i]]$y, pred[[1]])
        }
      })
    }
    mean_errs <- colMeans(cv.errs)
    min_ndims <- which(mean_errs == min(mean_errs))
    if(length(min_ndims) > 1){
      min_sd <- which.min(apply(cv.errs[,min_ndims], 2, sd))
      min_ndims <- min_ndims[min_sd]
    }
    out$dims <- out$dims[1:min_ndims]
    out$ProjectionMatrix <- out$ProjectionMatrix[,1:min_ndims]
    out$ProjectedData <- as.matrix(out$ProjectedData[,1:min_ndims])
    out$cv.errs <- mean_errs
  }

  ## Return
  if (!is.null(var_names)) {
    rownames(out$ProjectionMatrix) <- var_names
  }
  out
}
