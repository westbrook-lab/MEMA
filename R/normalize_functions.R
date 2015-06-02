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
#' @param baseECM A regular expressin for the name or names of the printed ECM(s) to be normalized against
#' @param baseGF A regular expression for the name or names of the soluble growth factors to be normalized against
#' @return A numeric vector of the normalized values
#'
#' @section Details: \code{normWellsWithinPlate} normalizes the value column of
#'   all MEPs by dividing the median value of the replicates of the MEP that
#'   is the pairing of baseECM  with baseGF.
#'   @export
normWellsWithinPlate <- function(DT, value, baseECM, baseGF) {
  if(!c("ShortName") %in% colnames(DT)) stop(paste("DT must contain a ShortNam column."))
  if(!c(value) %in% colnames(DT)) stop(paste("DT must contain a", value, "column."))
  if("Ligand" %in% colnames(DT)){
    valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ShortName)  & grepl(baseGF,DT$Ligand)),value, with=FALSE]), na.rm = TRUE)
  } else if (c("Growth.Factors") %in% colnames(DT)) {
    valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ShortName)  & grepl(baseGF,DT$Growth.Factors)),value, with=FALSE]), na.rm = TRUE)
  } else stop (paste("DT must contain a Growth.Factors or Ligand column."))
  normedValues <- DT[,value,with=FALSE]/valueMedian
  return(normedValues)
}

#' loessModel Create a median normalized loess model of an array
#'
#'@param data A dataframe with ArrayRow, ArrayColumn and signal intensity columns
#'@param value The column name of the signal intensity column
#'@param span The span value passed to loess. Values between 0 and 1 determine the
#'proportion of the population to be included in the loess neighborhood.
#'@return a vector of median normalized loess values of the signal
#'@export
loessModel <- function(data, value, span){
  dataModel <- loess(as.formula(paste0(value," ~ ArrayRow+ArrayColumn")), data,span=span)
  dataPredicted <- predict(dataModel)
  predictedMedian <- median(dataPredicted, na.rm = TRUE)
  dataNormed <- dataPredicted/predictedMedian
}
