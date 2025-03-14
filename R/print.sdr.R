#### Print for class "hldr" ####
print.sdr <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nsdrMethod:\n", paste0(x$method), "\n\n", sep = "")
  cat("\ndimselectMethod:\n",
      paste0(ifelse(is.null(x$dimselect), "none", x$dimselect)),
      "\n\n", sep = "")
  cat("\nProjectionMatrix:\n")
  print(x$ProjectionMatrix)
  invisible(x)
}
