#QA Functions for MEMAs

#' Randomize a vector of values
#' @param x a numeric vector
#'@return x with the values randomly ordered
#'@export
randomizePositions <- function(x){
  tmp <- x[sample(1:length(x),size = length(x), replace=FALSE)]
  return(tmp)
}


#' Calculate Coefficient of Variation (CV)
#'
#'Calculates the CV of a numeric vector as the standard deviation over the mean. All NA values are removed from the input vector.
#'@param x A numeric vector.
#'@return The CV of the non-na values in the x.
#'
#'  @export
CV <- function(x){
  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
}

# Empty function
#
cleanDT <- function (DT) {
  return(DT)
}
