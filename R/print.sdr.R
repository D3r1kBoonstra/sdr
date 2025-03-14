#### Print for class "hldr" ####
print.sdr <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nsdrMethod:\n", paste0(x$method), "\n\n", sep = "")
  cat("\ndim.orderMethod:\n",
      paste0(ifelse(is.null(x$dim.order.method), "none", x$dim.order.method)),
      "\n\n", sep = "")
  cat("\nProjectionMatrix:\n")
  print(x$ProjectionMatrix)
  invisible(x)
}
