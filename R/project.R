#### Projection Function ####
project <- function(x, ProjectionMatrix){
  if(!is.matrix(x)) x <- as.matrix(x)
  t(ProjectionMatrix %*% t(x))
}
