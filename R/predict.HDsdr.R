predict.HDsdr <- function(object, newdata){
  if(missing(newdata)){
    predict(MASS::qda(object$ProjectedData, object$slices))
  } else {
    newdata_center <- as.matrix(newdata) - object$center
    newdata_proj <- newdata_center %*% object$spca %*% object$ProjectionMatrix
    mod <- MASS::qda(object$ProjectedData, object$slices)
    predict(mod, newdata = newdata_proj)
  }
}
