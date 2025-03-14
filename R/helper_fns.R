S_W <- function(prior, data){
  Reduce(`+`,
         mapply(function(x, y){
           x * cov(y)
         }, prior, data, SIMPLIFY = FALSE)
  )
}
S_B <- function(prior, xbar){
  xbarbar <- Reduce(`+`,
                    mapply(function(x, y){x * y},
                           prior, xbar, SIMPLIFY = FALSE)
  )
  Reduce(`+`,
         mapply(function(x, y, z){
           (y - z) %*% t(y - z)
         }, prior, xbar, list(xbarbar = xbarbar), SIMPLIFY = FALSE)
  )
}

mat_power <-  function(x, power) {
  x <-  (x + t(x)) / 2
  with(
    eigen(x), vectors %*% diag(values^power) %*% t(vectors))

}

data_list_fn <- function(x, grouping){
  grouping <- droplevels(as.factor(grouping))
  classes <- unique(grouping)
  lapply(1:length(classes), function(i){
    x[grouping == classes[[i]], ]
  })
}

priors_fn <- function(x, grouping, data_list = NULL){
  if(is.null(data_list)) data_list <- data_list_f(x, grouping)
  n <- sapply(data_list, nrow)
  total <- sum(n)
  sapply(n, function(x){x / total})
}

lik <- function(d){
  p <- length(d)
  out <- vector(length = p)
  for(q in 1:(p-1)){
    d_list <- list(d[1:q], d[(q+1):p])
    mu <- lapply(d_list, mean)
    s2 <- lapply(d_list, function(x) ifelse(length(x) == 1, as.numeric(0), var(x)) )
    pooled_sd <- sqrt(((q-1)*s2[[1]] + (p-q-1)*s2[[2]] )/(p-2))

    out[[q]] <- sapply(seq_along(d_list), function(i){
      sum(dnorm(d_list[[i]], mu[[i]], pooled_sd, log = TRUE))
    }) |>
      sum()
  }
  out[[p]] <- sum(dnorm(d, mean(d), sqrt(var(d)), log = TRUE))
  ## Return
  list("lik" = out, "max" = which.max(out))
}

dot <- function(x, y) {
  stopifnot(length(x) == length(y))
  sum(Conj(y) * x)
}
norm <- function(x) sqrt(dot(x, x))
is.zero <- function(x) norm(x) < sqrt(.Machine$double.eps)

normalize <- function(x) {
  # if (is.zero(x)) {
  #   stop("Input `x` is numerically equivalent to the zero vector.", call. = FALSE)
  # } else {
  x / norm(x)
  # }
}
