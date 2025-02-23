#### Projection Function ####
project <- function(x, ProjectionMatrix){
  t(ProjectionMatrix %*% t(x))
}