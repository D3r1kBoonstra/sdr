#### hldr.fit function ####
sdr.fit <- function(x, grouping, method = "SDRS", dims = NULL, dimselect = NULL, ...){

  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.factor(grouping)) grouping <- as.factor(grouping)
  ## Removed unused factor levels
  grouping <- droplevels(grouping)
  if(is.null(dims)) dims <- 1:ncol(x)


  # Qda default code
  # if (is.null(dim(x)))
  #   stop("'x' is not a matrix")
  # x <- as.matrix(x)
  # if (any(!is.finite(x)))
  #   stop("infinite, NA or NaN values in 'x'")
  # n <- nrow(x)
  # p <- ncol(x)
  # if (n != length(grouping))
  #   stop("nrow(x) and length(grouping) are different")
  # g <- as.factor(grouping)
  # lev <- levels(g)
  # counts <- as.vector(table(g))
  # names(counts) <- lev
  # if (any(counts < p + 1))
  #   stop("some group is too small for 'qda'")
  # proportions <- counts/length(g)
  # ng <- length(proportions)
  # if (any(prior < 0) || round(sum(prior), 5) != 1)
  #   stop("invalid 'prior'")
  # if (length(prior) != ng)
  #   stop("'prior' is of incorrect length")


  if(method == "SDRS"){
    out <- sdrs(x, grouping, dims, ...)
  }
  if(method == "SDRS2"){
    out <- sdrs2(x, grouping, dims, ...)
  }
  if(method == "LD"){
    out <- LD(x, grouping, dims, ...)
  }
  if(method == "SIR"){
    out <- SIR(x, grouping, dims, ...)
  }
  if(method == "SAVE"){
    out <- SAVE(x, grouping, dims, ...)
  }

  if(!is.null(dimselect)) {
    dimselect_out <- dimselect(x, grouping, out$ProjectionMatrix,
                               method = dimselect, ...)
    out$ProjectionMatrix <- dimselect_out$ProjectionMatrix
  }

  class(out) <- c("sdr")

  out$ProjectedData <- as.matrix(project(x, out$ProjectionMatrix))
  out$grouping <- grouping
  out$dims <- dims

  # if(fitted.class){
  #   out$fitted.class <- predict(MASS::qda(out$ProjectedData, grouping), ...)
  # }
  if(!is.null(dimselect)){
    out$dims <- dimselect_out$dims
  }
  ## Return
  out
}
