#Functions to support MEMAs printed in 8 well plates


#'melt raw population dataset to combine wells
#'
#'\code{melt8Well}
#'
#'@param DT a data.table with columns of population intensity values in columns named by barcode, well, wavelength and data type(Raw/Net/Background)
#'@return A data.table in long format with columns of intensity readings for #' each data type and wavelength, barcode and well.
#'@section Usage: This function reformats Array Pro results files that
#' come from analyzing images from multiple channels, wells and plates.
#' The column names are encoded as follows:
#' Datatype is at the beginning of the column name as Raw, Net or
#' Background. Barcode follows the curly brace symbol. Wavelength follows
#'  the - symbol after the barcode. Well follows the - symbol after the
#'  wavelength.

#' @export
melt8Well<-function(DT){

  intNames<-grep("Background|Net|Raw",colnames(DT), value = TRUE)
  intNameType<-lapply(intNames,function(intName){
    return(strsplit(intName,split=" ")[[1]][1])
  })
  intNameWL<-lapply(intNames,function(intName){
    return(strsplit(intName,split="-")[[1]][2])
  })
  intNameWell<-lapply(intNames,function(intName){
    return(strsplit(intName,split="-")[[1]][3])
  })
  intNameBarcode<-lapply(intNames,function(intName){
    tmp<-strsplit(intName,split="[{]")[[1]][2]
    return(strsplit(tmp,split="-")[[1]][1])
  })
  newIntNames<-make.names(paste(intNameType,intNameWL,intNameWell,intNameBarcode))
  data.table::setnames(DT,colnames(DT[,intNames, with=FALSE]),newIntNames)
  meltDT<-reshape2::melt(DT,measure=newIntNames,variable="TypeWLWellBarcode",value="Intensity", variable.factor=FALSE)
  meltDT$Type<-unlist(lapply(meltDT$TypeWLWellBarcode,function(x){
    return(strsplit(x,split="[.]")[[1]][1])
  }))
  meltDT$WL<-unlist(lapply(meltDT$TypeWLWellBarcode,function(x){
    return(strsplit(x,split="[.]")[[1]][2])
  }))
  meltDT$Well<-unlist(lapply(meltDT$TypeWLWellBarcode,function(x){
    well<-strsplit(x,split="[.]")[[1]][3]
    wellChars<-strsplit(well,"")
    Well<-paste0(wellChars[[1]][1],"0",wellChars[[1]][2])
    return(Well)
  }))
  meltDT$Barcode<-unlist(lapply(meltDT$TypeWLWellBarcode,function(x){
    return(strsplit(x,split="[.]")[[1]][4])
  }))
  meltDT[,TypeWLWellBarcode:=NULL]
  meltDT$Type<-paste(meltDT$Type,meltDT$WL,sep=".")
  #Delete the Wavelength only column
  meltDT[,WL:=NULL]
  data.table::setnames(meltDT,colnames(meltDT),make.names(colnames(meltDT)))

  #wide format creating columns of data organized by type and wavelength
  DT<-reshape2::dcast(meltDT,...~Type,value.var="Intensity")
  DT<-data.table::data.table(DT,key=c("Well","Barcode"))
  return(DT)
}

#' Rotate the metadata 180 degrees in Array space
#' @export
rotateMetadata <- function(DT){
  DT$ArrayRow <- max(DT$ArrayRow)+1-DT$ArrayRow
  DT$ArrayColumn <- max(DT$ArrayColumn)+1-DT$ArrayColumn
  DT$Spot <- as.integer(max(DT$Spot)+1-DT$Spot)
  return(DT)
}
