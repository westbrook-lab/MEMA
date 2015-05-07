#Normalization functions for processing MEMAs

#' Normalize the proliferation ratio signal to the collagen 1 values
#' @param x a dataframe or datatable with columns names ProliferatioRatio
#' and ShortName. ShortName must include at least one entry of COL1 or COL I.
#' @return The input dataframe of datatable with a normedProliferation column that has the ProliferationRatio values divided by the median collagen
#' 1 proliferation value
#' @export
normProfToCol1 <- function(x){
  col1Median <- median(x$ProliferationRatio[x$ShortName %in% c("COL1", "COL I")],na.rm = TRUE)
  normedProliferation <- x$ProliferationRatio/col1Median
}


#' Normalize to a base MEP
#'
#' Normalizes one channel of values for all MEPs in a multi-well plate to one
#' base MEP.
#'
#' @param DT A \code{data.table} that includes a numeric value column to be
#'   normalized, a \code{ShortName} column that has the printed ECM names and a
#'   \code{Growth.Factors} column that has the growth factor names.
#' @param value The name of the column of values to be normalized
#' @param baseECM The name of the printed ECM to be normalized against
#' @param baseGF The name of the soluble growth factor to be normalized against
#' @return A numeric vector of the normalized values
#'
#' @section Details: \code{normWellsWithinPlate} normalizes the value column of
#'   all MEPs by subtracting the median value of the replicates of the MEP that
#'   is the pairing of baseECM  with baseGF.
normWellsWithinPlate <- function(DT, value, baseECM, baseGF) {
  if(!c("ShortName") %in% colnames(DT)) stop(paste("DT must contain a ShortNam column."))
  if(!c("Growth.Factors") %in% colnames(DT)) stop(paste("DT must contain a Growth.Factors column."))
  if(!c(value) %in% colnames(DT)) stop(paste("DT must contain a", value, "column."))
  valueMedian <- median(unlist(DT[(DT$ShortName == baseECM & DT$Growth.Factors == baseGF),value, with=FALSE]), na.rm = TRUE)
  normedValues <- DT[,value,with=FALSE]-valueMedian
  return(normedValues)
}

