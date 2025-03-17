#### Print for class "hldr" ####
print.sdr <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nsdrMethod:\n", paste0(x$method), "\n\n", sep = "")
  cat("\ndim.orderMethod:\n",
      paste0(ifelse(is.null(x$dim_order_method), "none", x$dim_order_method)),
      "\n\n", sep = "")
  if(!is.null(x$dim_order_method)){
    cat("\ndims:\n", paste0(x$dims), "\n\n", sep = "")
  }
  cat("\nProjectionMatrix:\n")
  print(x$ProjectionMatrix)
  invisible(x)
}
