# Helper function to create level 3 and 4 data
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))


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
  metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Endpoint|ECMp|MEP|Well_Ligand|ImageID|Barcode|^Well$|^PrintSpot$|^Spot$)", x=colnames(cDT), value=TRUE)
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
  l4Names<-grep("Loess$|RUV3|Norm|^Ligand|^ECMp|Barcode|Spot_PA_SpotCellCount$|Spot_PA_ReplicateCount$", x=names(l3),value=TRUE)
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
  
  l3Names <- grep("Barcode|Well|CellLine|Ligand|Endpoint488|Endpoint555|Endpoint647|EndpointDAPI|ECMp|MEP", colnames(l3), value=TRUE)
  #Merge back in the replicate metadata
  mDT <- l3[,l3Names,keyby="Ligand,ECMp,Barcode", with=FALSE]
  setkey(mDT,Ligand,ECMp,Barcode)
  l4DT <- mDT[l4DT, mult="first"]
  l4DT <- l4DTse[l4DT]
  
  return(l4DT)
}#End of createl4
