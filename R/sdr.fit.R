#### hldr.fit function ####
sdr.fit <- function(x, slices, method = "sdrs", dims = NULL, dim.order = NULL, ...){


  method <- tolower(method)
  if(method == "sdrs"){
    out <- sdrs(x, slices, dims, ...)
  }
  if(method == "sdrs2"){
    out <- sdrs2(x, slices, dims, ...)
  }
  if(method == "sdrs3"){
    out <- sdrs3(x, slices, dims, ...)
  }
  # if(method == "ld"){
  #   out <- LD(x, slices, dims, ...)
  # }
  if(method == "sir"){
    out <- sir(x, slices, dims, ...)
  }
  if(method == "save"){
    out <- save(x, slices, dims, ...)
  }
  if(method == "dr"){
    out <- dr(x, slices, dims, ...)
  }
  out$ProjectionMatrix <- as.matrix(out$ProjectionMatrix)
  class(out) <- c("sdr")

  out$ProjectedData <- x %*% out$ProjectionMatrix
  out$slices <- slices
  out$dims <- dims

  if(!is.null(dim.order)) {
    dim_order_out <- dim_order(object = out, method = dim.order, ...)
    out$ProjectionMatrix <- dim_order_out$ProjectionMatrix
    out$dims <- dim_order_out$dims
    out$ProjectedData <- dim_order_out$ProjectedData
    out$dim_criteria <- dim_order_out$dim_criteria
  }

  ## Return
  out
}
