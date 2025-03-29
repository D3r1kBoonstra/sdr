#### Print for class "sdr" ####
print.sdr <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nsdrMethod:\n", paste0(x$method), "\n\n", sep = "")
  if(!is.null(x$dim_order_method)){
    cat("\ndim.orderMethod:\n",
        paste0(ifelse(is.null(x$dim_order_method), "none", x$dim_order_method)),
        "\n\n", sep = "")
    cat("\ndims:\n", paste0(x$dims, sep = ","), "\n\n", sep = "")
    cat("\ndim_criteria:\n", x$dim_criteria, "\n\n", sep = " ")
  }
  cat("\nProjectionMatrix:\n")
  print(x$ProjectionMatrix)
  invisible(x)
}
