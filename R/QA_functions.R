#QA Functions for MEMAs

#' Randomize a vector of values
#' @param x a numeric vector
#'@return x with the values randomly ordered
#'@export
randomizePositions <- function(x){
  tmp <- x[sample(1:length(x),size = length(x), replace=FALSE)]
  return(tmp)
}
