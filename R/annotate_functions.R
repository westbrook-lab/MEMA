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

#'
#'@export
compressHA <- function(x){
  x <- gsub("(hyaluronic_acid_greater_than_500kDa)","HA>500kDa",x)
  x <- gsub("(hyaluronic_acid_less_than_500kDa)","HA<500kDa",x)
  x <- gsub("hyaluronicacid","HA",x)
  x <- gsub("lessthan","<",x)
  x <- gsub("greaterthan",">",x)
  return(x)
}

#'
#'@export
processan2omero <- function (fileNames) {
  rbindlist(lapply(fileNames, function(fn){
    #Process each file separately
    dt <- fread(fn,header = TRUE)
    
    #Rename to preprocessing pipeline variable names
    setnames(dt,"OSpot","Spot")
    setnames(dt,"PlateID","Barcode")
    setnames(dt,"395nm","EndpointDAPI")
    setnames(dt,"488nm","Endpoint488")
    setnames(dt,"555nm","Endpoint555")
    setnames(dt,"640nm","Endpoint647")
    setnames(dt,"750nm","Endpoint750")
    #Shorten and combine Annot names
    dt$CellLine <- gsub("_.*","",dt$CellLine)
    dt$ECM1 <- compressHA(dt$ECM1)
    dt$ECM2 <- compressHA(dt$ECM2)
    dt$ECM3 <- compressHA(dt$ECM3)
    #Chain ECM proteins if the second one is not COL1
    dt$ECMp <-paste0(gsub("_.*","",dt$ECM1),"_",gsub("_.*","",dt$ECM2),"_",gsub("_.*","",dt$ECM3)) %>%
      gsub("_NA","",.) %>%
      gsub("_COL1|_$","",.)
    #Chain ligands
    dt$Ligand <-paste0(gsub("_.*","",dt$Ligand1),"_",gsub("_.*","",dt$Ligand2)) %>%
      gsub("_NA","",.)
    dt$MEP <- paste0(dt$ECMp,"_",dt$Ligand)
    dt$Drug <- gsub("_.*","",dt$Drug1)
    dt$MEP_Drug <-paste0(dt$MEP,"_",dt$Drug)
    dt$EndpointDAPI <-gsub("_.*","",dt$EndpointDAPI)
    dt$Endpoint488 <-gsub("_.*","",dt$Endpoint488)
    dt$Endpoint555 <-gsub("_.*","",dt$Endpoint555)
    dt$Endpoint647 <-gsub("_.*","",dt$Endpoint647)
    #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
    dt$PrintSpot <- dt$Spot
    nrArrayRows <- max(dt$ArrayRow)
    nrArrayColumns <- max(dt$ArrayColumn)
    dt$PrintSpot[grepl("B", dt$Well)] <- (nrArrayRows*nrArrayColumns+1)-dt$PrintSpot[grepl("B", dt$Well)]
    return(dt)
    #  }, mc.cores=max(4, detectCores())))
  }))
  
}

