make_slices <- function(y, nslices = 5, ...) {
  if(is.character(y) || is.factor(y)){
    warning(paste0("y is categorical.", " 'nslices = ", nslices, "' is being ignored. Returning slices by defined categories."))
    if(!is.factor(y)) y <- as.factor(y)
    out <- droplevels(y)
  } else if(is.numeric(y)){
    n <- length(y)
    if (nslices > n) stop(paste0("nslices cannot be greater than the length of y."), "'nslices = ", nslices, "' and length(y) = ", n, ".")
    m <- floor(n / nslices)

    y_ord <- sort(y)
    slice_points <- y_ord[seq(m + 1, m * (nslices - 1) + 1, by = m)]

    out <- as.factor(findInterval(y, slice_points) + 1)
  } else {
    stop(paste0("y must be numeric, character, or factor. typeof(y) = ", typeof(y), "."))
  }
  out
}
