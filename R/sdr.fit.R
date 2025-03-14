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

  if(!is.null(dim.order)) {
    dim_order_out <- dim_order(x, slices, out$ProjectionMatrix,
                               method = dim.order, ...)
    out$ProjectionMatrix <- dim_order_out$ProjectionMatrix
  }

  class(out) <- c("sdr")

  out$ProjectedData <- x %*% out$ProjectionMatrix
  out$slices <- slices
  out$dims <- dims

  # if(fitted.class){
  #   out$fitted.class <- predict(MASS::qda(out$ProjectedData, slices), ...)
  # }
  if(!is.null(dim.order)){
    out$dims <- dim_order_out$dims
  }
  ## Return
  out
}
