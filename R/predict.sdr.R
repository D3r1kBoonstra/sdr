#### Predict function for class "hldr" ####
predict.sdr <- function(object, newdata = NULL, ndims = NULL, type = "model", model = "guess", ...){
  if(is.null(ndims)) ndims <- ncol(as.matrix(object$ProjectionMatrix))
  dims <- 1:ndims
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
      if(model == "guess"){
        if(object$ytype == "categorical"){
          model <- MASS::qda
        } else {
          model <- lm
        }
      } else if(model == "qda"){
        model <- MASS::qda
      } else if(model == "lda") {
        model <- MASS::lda
      } else if(model == "rf"){
        model <- randomForest::randomForest
      } else if(model == "lm"){
        model <- lm
      }
    }

    if(is.null(newdata)){
      df <- data.frame("y" = object$y, as.matrix(object$ProjectedData[,dims]))
      out <- predict(model(y ~., data = df,  ...))
    } else if(!is.null(newdata)){
      train_proj <- as.matrix(object$ProjectedData[,dims])
      colnames(train_proj) <- paste0("dim", dims)
      train_df <- data.frame("y" = object$y, train_proj)
      test_proj <- as.data.frame(as.matrix(newdata) %*% ProjectionMatrix)
      colnames(test_proj) <- colnames(train_proj)
      train_mod <- model(y ~., data = train_df, ...)
      out <- predict(train_mod, newdata = test_proj)
    }
  }
  if(!is.list(out)) out <- list(out)
  out$ProjectionMatrix <- ProjectionMatrix
  out

}
