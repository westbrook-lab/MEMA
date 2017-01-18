# Helper functions to create level 3 and 4 data

#' Calculate standard error of the mean
#' 
#' Omit na values if present
#' @param x a numeric vector
#' @return the standard erro of the mean for x as a numeric value
#' @export
se <- function(x) sd(x)/sqrt(sum(!is.na(x)))

#' Summarize cell level data to the spot level
#' 
#' Median summarize the cell level normalized values and the most biologically 
#' interpretable raw data at each spot, calculate the standard errors and add
#' SE columns for all median summarized data
#' @param cDT The datatable of cell level data to be summarized
#' @param lthresh The threshold used in the loess model to define low cell count
#'  regions
#'  @return A datatable of spot-level, median summarized data with standard error values and 
#'  metadata
#' @export
createl3 <- function(cDT, lthresh = lthresh, seNames=NULL){
  #Summarize cell data to medians of the spot parameters
  parameterNames<-grep(pattern="(Children|_CP_|_PA_)",x=names(cDT),value=TRUE)
  
  #Remove spot-normalized or summarized parameters
  parameterNames <- grep("SpotNorm|Wedge|Sparse|OuterCell|Center|^Nuclei_PA_Gated_EduPositive$",parameterNames,value=TRUE,invert=TRUE)
  
  slDT<-cDT[,lapply(.SD, numericMedian), by="Barcode,Well,Spot", .SDcols=parameterNames]
  #Use seNames to select the parameters that get SE values
  if(!is.null(seNames)){
    seNamesPattern<-paste(seNames,collapse="|")
    seNames <- grep(seNamesPattern,parameterNames,value=TRUE)
    slDTse <- cDT[,lapply(.SD,se), by="Barcode,Well,Spot", .SDcols=seNames]
  } else{
    slDTse <- cDT[,lapply(.SD,se), by="Barcode,Well,Spot"]
  }
  
  #Add _SE to the standard error column names
  setnames(slDTse, grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the spot and well metadata
  metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Endpoint|ECMp|MEP|Well_Ligand|ECM|ImageID|Barcode|^Well$|^PrintSpot$|^Spot$|Pin|Lx)", x=colnames(cDT), value=TRUE)
  setkey(cDT,Barcode, Well,Spot)
  mDT <- cDT[,metadataNames,keyby="Barcode,Well,Spot", with=FALSE]
  slDT <- mDT[slDT, mult="first"]
  #Merge in the standard err values
  setkey(slDTse, Barcode, Well, Spot)
  slDT <- slDTse[slDT]
  #Add a count of replicates
  slDT <- slDT[,Spot_PA_ReplicateCount := .N,by="Ligand,ECMp"]
  
  #Add the loess model of the SpotCellCount on a per well basis
  slDT <- slDT[,Spot_PA_LoessSCC := loessModel(.SD, value="Spot_PA_SpotCellCount", span=.5), by="Barcode,Well"]
  
  #Add well level QA Scores to spot level data
  slDT <- slDT[,QAScore := calcQAScore(.SD, threshold=lthresh, maxNrSpot = max(cDT$ArrayRow)*max(cDT$ArrayColumn),value="Spot_PA_LoessSCC"),by="Barcode,Well"]
  return(slDT)
}



numericMedianUniqueMetadata<-function(x){
  if(is.numeric(x)){
    as.numeric(median(x))
  } else {
    if(!length(unique(x))==1){
      return(NA)
    } else{
      unique(x)
    }
  }
}

#' @export
summarizeFBS <- function(dt){
  #Summarize all by ligand and ECMp
  #This will only change the FBS data
  dtFBS <- dt[grepl("FBS",dt$Ligand),]
  #Remove anything after FBS
  dtFBS$Ligand <- gsub("FBS.*","FBS", dtFBS$Ligand)
  #Create the new MEP name
  dtFBS$MEP <- paste(dtFBS$ECMp,dtFBS$Ligand,sep="_")
  if("StainingSet" %in% colnames(dtFBS)){
    #Find the medians or the unique metadata values
    dtFBSMedians <- dtFBS[, lapply(.SD, numericMedianUniqueMetadata), keyby = "Ligand,ECMp,StainingSet"]
  } else {
    #Find the medians or the unique metadata values
    dtFBSMedians <- dtFBS[, lapply(.SD, numericMedianUniqueMetadata), keyby = "Ligand,ECMp"]
  }
  #If the FBS is from multiple plates its barcode has been reset to NA.
  #replace this with the word Multiple
  dtFBSMedians$Barcode[is.na(dtFBSMedians$Barcode)] <- "Multiple"
  #Delete the FBS wells from the original dt
  dt <- dt[!grepl("FBS",dt$Ligand),]
  #Bind in the summarized FBS values
  dt <- rbind(dt,dtFBSMedians)
  return(dt)
}

#' Summarize spot level data to the MEP level
#' 
#' Median summarize the spot level normalized values the most biologically 
#' interpretable raw data at each spot, calculate the standard errors and add
#' SE columns for all median summarized data
#' @param l3 The datatable of spot level data to be summarized
#' @return A datatable of MEP level, median summarized data with standard error values and 
#'  metadata
#' @export
createl4 <- function(l3, seNames=NULL){
  #Add a count of replicates
  l3 <- l3[,Spot_PA_ReplicateCount := .N,by="Ligand,ECMp"]
  l4Names<-grep("Loess$|RUV|Norm|^Ligand$|^ECMp|Barcode|Spot_PA_SpotCellCount$|Spot_PA_ReplicateCount$", x=names(l3),value=TRUE)
  #remove the _SE values
  l4Names <- grep("_SE|NormMethod|AnnotID",l4Names, value = TRUE, invert = TRUE)
  l4Keep<-l3[,l4Names,with=FALSE]
  l4DT<-l4Keep[,lapply(.SD,numericMedian),keyby="Ligand,ECMp,Barcode"]
  #Use seNames to select the parameters that get SE values
  if(!is.null(seNames)){
    seNamesPattern<-paste(seNames,collapse="|")
    seNames <- grep(seNamesPattern,l4Names,value=TRUE)
    l4DTse <- l4Keep[,lapply(.SD,se),keyby="Ligand,ECMp,Barcode", .SDcols=seNames]
  } else{
    l4DTse <- l4Keep[,lapply(.SD,se),keyby="Ligand,ECMp,Barcode"]
  }
  
  #Add _SE to the standard error column names
  setnames(l4DTse, grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE),"_SE"))
  
  l3Names <- grep("Barcode|Well|CellLine|Ligand|ECM|Endpoint488|Endpoint555|Endpoint647|EndpointDAPI|ECMp|MEP|Lx|PinDiameter", colnames(l3), value=TRUE)
  #Merge back in the replicate metadata
  mDT <- l3[,l3Names,keyby="Ligand,ECMp,Barcode", with=FALSE]
  setkey(mDT,Ligand,ECMp,Barcode)
  l4DT <- mDT[l4DT, mult="first"]
  l4DT <- l4DTse[l4DT]
  l4DT <- summarizeFBS(l4DT)
  return(l4DT)
}#End of createl4

