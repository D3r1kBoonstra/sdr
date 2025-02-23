#### Predict function for class "hldr" ####
predict.sdr <- function(object, newdata = NULL, dims, grouping){
  if(missing(newdata)){
    # class <- model.frame(object)[,1]
    predict(MASS::qda(object$ProjectedData, grouping = object$grouping))
  } else if(!is.null(newdata)){
    if(!missing(dims))
    # dims <- 1:nrow(object$ProjectionMatrix)
    ProjectionMatrix <- object$ProjectionMatrix[dims,]
    train_proj <- as.matrix(object$ProjectedData[,dims])
    test_proj <- as.matrix(project(newdata, ProjectionMatrix))
    train_qda <- MASS::qda(train_proj, grouping = object$grouping)
    out <- predict(train_qda, newdata = test_proj)
    out$ProjectionMatrix <- ProjectionMatrix
    out
  }
}
