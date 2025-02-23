predict.HDsdr <- function(object, newdata){
  if(missing(newdata)){
    predict(MASS::qda(object$ProjectedData, object$grouping))
  } else {
    newdata_center <- as.matrix(newdata) - object$center
    newdata_proj <- project(newdata_center %*% object$spca, object$ProjectionMatrix)
    mod <- MASS::qda(object$ProjectedData, object$grouping)
    predict(mod, newdata = newdata_proj)
  }
}
