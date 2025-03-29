#### hldr.fit function ####
sdr.fit <- function(method, ...){

  method <- tolower(method)
  class(method) <- method
  UseMethod("sdr.fit", method)


  # out$ProjectionMatrix <- as.matrix(out$ProjectionMatrix[,dims])
  # out$ProjectedData <- x %*% out$ProjectionMatrix
  # out$slices <- slices
  # class(out) <- c("sdr")

  # if(!is.null(dim.order)) {
  #   dim_order_out <- dim_order(object = out, method = dim.order, ndims = ndims, ...)
  #   out$ProjectionMatrix <- dim_order_out$ProjectionMatrix
  #   out$dims <- dim_order_out$dims
  #   out$eigvalues <- out$eigvalues[dims]
  #   out$ProjectedData <- dim_order_out$ProjectedData
  #   out$dim_criteria <- dim_order_out$dim_criteria
  # } else{
  #   dims <- dims[1:ndims]
  #   out$eigvalues <- out$eigvalues[dims]
  #   out$dims <- dims
  #   out$ProjectionMatrix <- as.matrix(out$ProjectionMatrix[,dims])
  #   out$ProjectedData <- out$ProjectedData[,dims]
  # }

  # out$ProjectionMatrix <- as.matrix(out$ProjectionMatrix[,dims])

  ## Return
  # out
}
