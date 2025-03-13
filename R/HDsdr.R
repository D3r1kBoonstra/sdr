HDsdr <- function(x, grouping, method = "SDRS2", dims = NULL, dim.order = NULL, cv.folds = NULL, orth = TRUE, spca.prop.var = TRUE, center = TRUE, sumabsv = NULL, npca = NULL, ...){

  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.factor(grouping)) grouping <- as.factor(grouping)
  grouping <- droplevels(grouping)
  if(is.null(sumabsv)){
    p <- ncol(x)
    max_sumv <- floor(sqrt(p))

    if(p < 5000L){
      if(max_sumv <= 20){
        sum_seqv <- seq(1, max_sumv, by = 1)
      } else {
        sum_seqv <- c(seq(0, max_sumv, by = 5), max_sumv)
        sum_seqv[1] <- 1
      }
    } else {
      sum_seqv <- seq(1, max_sumv, length.out = 10)
    }
    sumabsv <- PMA::SPC.cv(x, sumabsvs = sum_seqv, orth = orth, center = center, trace = FALSE)$bestsumabsv
  }
  if(is.null(npca)){
    classes <- unique(grouping)
    ni <- sapply(1:length(classes), function(i){
      nrow(x[grouping == classes[[i]], ])
    })
    npca <- min(ni) - 1
    if(!is.null(cv.folds)){
      npca <- floor(min(ni)*(1 - cv.folds/100) - 1)
    }
  }
  # if(center){
  #   xbar <- colMeans(x)
  #   scale <- apply(x, 2, sd)
  #   x <- x - xbar
  # }

  x_spca <- PMA::SPC(x, sumabsv = sumabsv, orth = orth, K  = npca, trace = FALSE, compute.pve = spca.prop.var, center = center)
  center <- ifelse(is.null(x_spca$meanx), 0, x_spca$meanx)
  out <- sdr((x - center) %*% x_spca$v, grouping, method, dims, dim.order, cv.folds, ... = ...)


  class(out) <- c("HDsdr")
  out$call <- match.call()
  out$SDRmethod <- method
  out$spca <- x_spca$v
  out$sumabsv <- sumabsv
  out$npca <- npca
  out$center <- center
  # out$scale <- scale
  if(spca.prop.var) out$spca.prop.var <- x_spca$prop.var.explained
  # # out$dimselect <- dimselect
  # # out$dims <- if(!is.null(dimselect)) out$dims
  # out$grouping <- grouping
  ## Return
  out
}
