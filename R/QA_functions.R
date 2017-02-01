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

#' Calculate a QA score for MEMAs
#'
#' Return a QA score based on the number of spots missing or below a threshold
#' @param DT A datatable with MEMA data
#' @param threshold A numeric vector of length 1
#' @param maxNrSpot A numeric vector of length 1 that is the maximum number of spots printed on the MEMA
#' @param value A character vector of length 1 that is the column name of the values that determine the QA score
#' @return A single numeric value that is the proportion of spots below the threshold
#'
#' @export
calcQAScore <- function(DT, threshold, maxNrSpot=700, value){
  QAScore <- (nrow(DT)-sum(DT[,value,with=FALSE] < threshold))/maxNrSpot
  return (QAScore)
}

#'
#' @export
RZScore <- function(x){
  xMedian <- median(x, na.rm=TRUE)
  xMad <-mad(x, na.rm=TRUE)
  if(xMad == 0){ zscores <- NA
  } else zscores <- (x-xMedian)/xMad
  return(zscores)
}
