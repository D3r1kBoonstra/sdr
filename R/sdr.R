sdr <- function(x, y, method = "sdrs", ytype = "guess", nslices = 5, ndims = NULL, dims = NULL, dim.order = NULL, cv.folds = NULL, ...){

  if(!is.matrix(x)) x <- as.matrix(x)

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
    if (is.factor(y) || is.character(y) || (length(unique(y)) / length(y) > 0.05)) {
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


  out <- sdr.fit(method, x = x, y = y, ytype = ytype, dims = dims, ...)
  class(out) <- c("sdr")
  out$call <- match.call()
  out$x <- x
  out$y <- y
  out$ytype <- ytype
  out$method <- method
  out$dims <- dims
  # out$ProjectionMatrix <- as.matrix(out$ProjectionMatrix[,dims])
  # out$ProjectedData <- x %*% out$ProjectionMatrix


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

  # if(!is.null(cv.folds)){
  #   df <- cbind("class" = slices, as.data.frame(x))
  #   splits <- rsample::vfold_cv(df, strata = class, v = cv.folds)$splits
  #   train <- lapply(splits, rsample::training)
  #   test <- lapply(splits, rsample::testing)
  #   cv.errs <- matrix(nrow = cv.folds, ncol = length(out$dims))
  #   for(i in 1:cv.folds){
  #     fit <- sdr.fit(x = as.matrix(train[[i]][,-1]), slices = train[[i]]$class, method = method, dims = out$dims, dim.order = NULL, lam = out$lam, ...)
  #     cv.errs[i,] <- sapply(1:length(out$dims), function(j){
  #       pred <- predict(fit, newdata = test[[i]][,-1], dims = out$dims[1:j], ...)$class
  #       mean(pred != test[[i]]$class)
  #     })
  #   }
  #   mean_errs <- colMeans(cv.errs)
  #   min_dims <- which(mean_errs == min(mean_errs))
  #   if(length(min_dims) > 1){
  #     min_sd <- which.min(apply(cv.errs[,min_dims], 2, sd))
  #     min_dims <- min_dims[min_sd]
  #   }
  #   out$dims <- out$dims[1:min_dims]
  #   out$ProjectionMatrix <- out$ProjectionMatrix[,1:min_dims]
  #   out$ProjectedData <- as.matrix(out$ProjectedData[,1:min_dims])
  #   out$cv.errs <- mean_errs
  # }

  ## Return
  if (!is.null(var_names)) {
    rownames(out$ProjectionMatrix) <- var_names
  }
  out
}
