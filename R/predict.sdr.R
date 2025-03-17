#### Predict function for class "hldr" ####
predict.sdr <- function(object, newdata = NULL, dims = NULL, type = "model", model = "qda", ...){
  if(is.null(dims)) dims <- 1:ncol(object$ProjectionMatrix)
  ProjectionMatrix <- as.matrix(object$ProjectionMatrix[,dims])

  if(type == "project") {
    out <- vector(mode = "list")
    if(is.null(newdata)){
      out$ProjectedData <- object$ProjectedData[,dims]
    } else {
      out$ProjectedData <- as.matrix(newdata) %*% ProjectionMatrix
    }
  }
  if(type == "model"){

    if(!is.function(model)){
      model <- tolower(model)
      if(model == "qda"){
        model <- MASS::qda
      } else if(model == "lda") {
        model <- MASS::lda
      } else if(model == "rf"){
        model <- randomForest::randomForest
      }
    }


    if(is.null(newdata)){
      out <- predict(model(as.matrix(object$ProjectedData[,dims]), object$slices, ...))
    } else if(!is.null(newdata)){
      train_proj <- as.matrix(object$ProjectedData[,dims])
      test_proj <- as.matrix(newdata) %*% ProjectionMatrix
      train_mod <- model(train_proj, object$slices, ...)
      out <- predict(train_mod, newdata = test_proj)
    }
  }
  if(!is.list(out)) out <- list(out)
  out$ProjectionMatrix <- ProjectionMatrix
  out

}
