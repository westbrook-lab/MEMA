#Functions to support MEMAs printed in 8 well plates


#' Rotate the metadata 180 degrees in Array space
#'
#'@param DT A data.table of metadata with Spot, ArrayRow and ArrayColumn columns.
#'@return The same data.table rotated 180 degrees in array space. The ArrayRow, Arraycolumn and Spot values are updated.
#'
#' @export
rotateMetadata <- function(DT){
  DT$ArrayRow <- max(DT$ArrayRow)+1-DT$ArrayRow
  DT$ArrayColumn <- max(DT$ArrayColumn)+1-DT$ArrayColumn
  DT$Spot <- as.integer(max(DT$Spot)+1-DT$Spot)
  return(DT)
}

#' Read in and parse an Aushon XML log file
#'
#' @param logFile An Aushon logfile
#' @return A datatable keyed by Row and Column with Depositions and
#'  PrintOrder columns.
#'
#' @export
readLogData<-function(logFile){
  #browser()
  data<-XML::xmlParse(logFile)
  dataList<-XML::xmlToList(data)
  #Only keep the sample attributes
  dataList<-dataList[names(dataList)=="Sample"]
  #Bind the XML data into a data table
  data<-data.table::rbindlist(dataList)
  #Create Row and Column data by shifting the values by 1
  data$Row<-as.integer(data$SpotRow)+1
  data$Column<-as.integer(data$SpotCol)+1
  #Convert deposition to an integer
  data$Depositions<-as.integer(data$Depositions)
  #Remove the 0 deposition entries
  data<-data[data$Depositions!=0,]
  #Create a print order column
  data$PrintOrder<-1:nrow(data)
  data.table::setkey(data,"PrintOrder")
  #Remove unneeded columns
  data <- data[,c("Row","Column","PrintOrder","Depositions"), with=FALSE]
  #Rotate by 90 degrees CW to match gal file orientation
  tmp <- data$Row
  data$Row <- data$Column
  data$Column <- 1+max(tmp)-tmp
  DT <- data.table::data.table(data,key="Row,Column")
  return(DT)
}

#' Convert column names in a data.table
#'
#' @param DT A data.table
#'
#' @return DT The same data.table with duplicated columns, invalid column name characters and trailing spaces removed.
#'
#' @export
convertColumnNames <- function (DT) {
  #Delete any duplicate names keeping the first instance
  DT <- DT[, unique(colnames(DT)), with = FALSE]
  #Replace invalid characters with a '.'
  data.table::setnames(DT, colnames(DT), make.names(colnames(DT)))
  #Remove all '.'s
  data.table::setnames(DT, colnames(DT), gsub("[.]", "", colnames(DT)))
}

#' Return the median of a vector as a numeric value
#'
#'\code{numericMedian} is a helper function for use within data.table that ensure the all medians are returned as numeric instead of numeric or integer values.
#' @param x integer or double vector
#'
#' @return The median of x as a numeric value
#'
#'@export
numericMedian <- function(x) as.numeric(median(x))
