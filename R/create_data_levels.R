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
#'  @export
createl3 <- function(cDT, lthresh = lthresh){
  #Summarize cell data to medians of the spot parameters
  parameterNames<-grep(pattern="(Children|_CP_|_PA_|Barcode|^Spot$|^Well$)",x=names(cDT),value=TRUE)
  
  #Remove any spot-normalized and cell level parameters
  parameterNames <- grep("SpotNorm|^Nuclei_PA_Gated_EduPositive$|^Nuclei_PA_Gated_EduPositive_RZSNorm$",parameterNames,value=TRUE,invert=TRUE)
  
  #Remove any raw parameters
  parameterNames <- grep("Barcode|^Spot$|^Well$|Norm|Nuclei_CP_Intensity_MedianIntensity_Dapi$|Cytoplasm_CP_Intensity_MedianIntensity_Actin$|Cytoplasm_CP_Intensity_MedianIntensity_CellMask$|Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker$|Nuclei_CP_Intensity_MedianIntensity_H3$|Nuclei_CP_Intensity_MedianIntensity_Fibrillarin$|Nuclei_CP_Intensity_MedianIntensity_Edu$|Cytoplasm_CP_Intensity_MedianIntensity_KRT5$|Cytoplasm_CP_Intensity_MedianIntensity_KRT19$|Spot_PA_SpotCellCount$", parameterNames, value = TRUE)
  
  cDTParameters<-cDT[,parameterNames,with=FALSE]
  
  slDT<-cDTParameters[,lapply(.SD,numericMedian),keyby="Barcode,Well,Spot"]
  slDTse <- cDTParameters[,lapply(.SD,se),keyby="Barcode,Well,Spot"]
  
  #Add _SE to the standard error column names
  setnames(slDTse, grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the spot and well metadata
  #TODO: Convert the logic to not name the metadata
  metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Endpoint|ECMp|MEP|Barcode|^Well$|^Spot$)", x=colnames(cDT), value=TRUE)
  setkey(cDT,Barcode, Well,Spot)
  mDT <- cDT[,metadataNames,keyby="Barcode,Well,Spot", with=FALSE]
  slDT <- mDT[slDT, mult="first"]
  #Merge in the standard err values
  slDT <- slDTse[slDT]
  #Add a count of replicates
  slDT <- slDT[,Spot_PA_ReplicateCount := .N,by="LigandAnnotID,ECMpAnnotID"]
  
  #Add the loess model of the SpotCellCount on a per well basis
  slDT <- slDT[,Spot_PA_LoessSCC := loessModel(.SD, value="Spot_PA_SpotCellCount", span=.5), by="Barcode,Well"]
  
  #Add well level QA Scores to spot level data
  slDT <- slDT[,QAScore := calcQAScore(.SD, threshold=lthresh, maxNrSpot = max(cDT$ArrayRow)*max(cDT$ArrayColumn),value="Spot_PA_LoessSCC"),by="Barcode,Well"]
}#End of create l3


#' Summarize spot level data to the MEP level
#' 
#' Median summarize the spot level normalized values the most biologically 
#' interpretable raw data at each spot, calculate the standard errors and add
#' SE columns for all median summarized data
#' @param l3 The datatable of spot level data to be summarized
#' @return A datatable of MEP level, median summarized data with standard error values and 
#'  metadata
#' @export
createl4 <- function(l3){
  #Add a count of replicates
  l3 <- l3[,Spot_PA_ReplicateCount := .N,by="LigandAnnotID,ECMpAnnotID"]
  l4Names<-grep("Norm|LigandAnnotID|ECMpAnnotID|Barcode|Spot_PA_SpotCellCount$|Spot_PA_ReplicateCount$", x=names(l3),value=TRUE)
  #remove the _SE values
  l4Names <- grep("_SE",l4Names, value = TRUE, invert = TRUE)
  l4Keep<-l3[,l4Names,with=FALSE]
  l4DT<-l4Keep[,lapply(.SD,numericMedian),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
  
  l4DTse <- l4Keep[,lapply(.SD,se),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
  #Add _SE to the standard error column names
  setnames(l4DTse, grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the replicate metadata
  mDT <- l3[,list(Well,CellLine,Ligand,Endpoint488,Endpoint555,Endpoint647,EndpointDAPI,ECMp,MEP),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
  l4DT <- mDT[l4DT, mult="first"]
  l4DT <- l4DTse[l4DT]
  
  return(l4DT)
}#End of createl4
