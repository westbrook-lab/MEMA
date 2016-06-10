#'Apply RUV3 and Loess Normalization to the common signals in a dataset
#' @export
normRUV3LoessResidualsCS <- function(dt, k){
  setkey(dt,Barcode,Well,Ligand,ECMp)
  signalNames <- grep("_CP_|_PA_",colnames(dt), value=TRUE)
  
  #Add residuals from subtracting the biological medians from each value
  residuals <- dt[,lapply(.SD,calcResidual), by="Barcode,Well,Ligand,ECMp", .SDcols=signalNames]
  #Add within array location metadata
  residuals$Spot <- as.integer(dt$Spot)
  residuals$PrintSpot <- as.integer(dt$PrintSpot)
  residuals$ArrayRow <- dt$ArrayRow
  residuals$ArrayColumn <- dt$ArrayColumn
  #Create a signal type
  dt$SignalType <- "Signal"
  residuals$SignalType <- "Residual"
  srDT <- rbind(dt,residuals)
  
  #Add to carry metadata into matrices
  srDT$BWL <- paste(srDT$Barcode, srDT$Well, srDT$Ligand, sep="_") 
  #srDT$SERC <- paste(srDT$Spot,srDT$ECMp, srDT$ArrayRow, srDT$ArrayColumn, sep="_")
  #srDT$SRC <- paste(srDT$Spot, srDT$ArrayRow, srDT$ArrayColumn, sep="_")
  
  #Set up the M Matrix to denote replicates
  nrControlWells <- sum(grepl("FBS",unique(srDT$BWL[srDT$SignalType=="Signal"])))
  nrLigandWells <- length(unique(srDT$Ligand[srDT$SignalType=="Signal"]))-nrControlWells
  M <-matrix(0, nrow = length(unique(srDT$BWL[srDT$SignalType=="Signal"])), ncol = nrLigandWells+1)
  rownames(M) <- unique(srDT$BWL[srDT$SignalType=="Signal"])
  #Indicate the control wells in the last column
  Mc <- M[grepl("FBS",rownames(M)),]
  Mc[,ncol(Mc)] <-1L
  #Subset to the ligand wells and mark as non-replicate
  Ml <- M[!grepl("FBS",rownames(M)),]
  #Name the rows by their ligand names, removing barcode and well
  MlBWLNames <- rownames(Ml)
  rownames(Ml) <- sub(".*?_","",sub(".*?_","",rownames(Ml)))
  #Name the columns a unique ligand + FBS
  colnames(Ml)<- c(unique(rownames(Ml)),"FBS")
  MlList <- lapply(colnames(Ml),function(cn,Ml){
    #Take the ligand M matrix and put a 1 where the rowname==colname
    Ml[rownames(Ml)==cn,colnames(Ml)==cn] <- 1
    return(Ml[,colnames(Ml)==cn])
  },Ml=Ml)
  Ml <- do.call(cbind, MlList)
  #Restore the row names ot the full BWL values
  rownames(Ml)<- MlBWLNames
  #Add the replicate wells and restore the row order
  M <- rbind(Mc,Ml)
  M <- M[order(rownames(M)),]
  
  srmList <- lapply(signalNames, function(signalName, dt){
    srm <- MEMA:::signalResidualMatrix(dt[,.SD, .SDcols=c("BWL", "PrintSpot", "SignalType", signalName)])
    return(srm)
  },dt=srDT)
  
  names(srmList) <- signalNames
  
  srmRUV3List <- mclapply(names(srmList), function(srmName, srmList, M, k){
    Y <- srmList[[srmName]]
    #Hardcode in identification of residuals as the controls
    resStart <- ncol(Y)/2+1
    cIdx=resStart:ncol(Y)
    nY <- RUVIIIArrayWithResiduals(k, Y, M, cIdx, srmName) #Normalize the spot level data
    nY$k <- k
    nY$SignalName <- paste0(srmName,"RUV3")
    setnames(nY,srmName,paste0(srmName,"RUV3"))
    return(nY)
  }, srmList=srmList, M=M, k=k, mc.cores=detectCores())
  
  #Reannotate with ECMp and MEP
  ECMpDT <- unique(srDT[,list(Well,PrintSpot,Spot,ECMp,ArrayRow,ArrayColumn)])
  
  srmERUV3List <- mclapply(srmRUV3List, function(dt,ECMpDT){
    setkey(dt,Well,PrintSpot)
    setkey(ECMpDT,Well,PrintSpot)
    dtECMpDT <- merge(dt,ECMpDT)
    dtECMpDT$MEP <- paste(dtECMpDT$ECMp,dtECMpDT$Ligand,sep="_")
    return(dtECMpDT)
  },ECMpDT=ECMpDT, mc.cores=detectCores())
  
  
  #Add Loess normalized values for each signal
  RUV3LoessList <- mclapply(srmERUV3List, function(dt){
    dtRUV3Loess <- loessNormArray(dt)
  },mc.cores=detectCores())
  
  #Combine the normalized signal into one data.table
  #with one set of metadata
  signalDT <- do.call(cbind,lapply(RUV3LoessList, function(dt){
    sdt <- dt[,grep("_CP_|_PA_",colnames(dt)), with=FALSE]
  }))
  
  signalMetadataDT <- cbind(RUV3LoessList[[1]][,grep("_CP_|_PA_",colnames(RUV3LoessList[[1]]), invert=TRUE), with=FALSE], signalDT)
  signalMetadataDT <- signalMetadataDT[,SignalName := NULL]
  signalMetadataDT <- signalMetadataDT[,mel := NULL]
  signalMetadataDT <- signalMetadataDT[,Residual := NULL]
  
  #Merge the raw signal back into the dataset
  signalMetadataDT <- merge(signalMetadataDT,dt,by=c("Barcode","Well","Spot","PrintSpot","ArrayRow","ArrayColumn","ECMp","Ligand"))
  
  return(signalMetadataDT)
}

#'Apply RUV3 and Loess Normalization to the common signals in a dataset
#' @export
normRUV3LoessResidualsCSx2 <- function(dt, k){
  #dt <- l3n[,colnames(l3n)[1:13], with=FALSE]
  setkey(dt,Barcode,Well,Ligand,ECMp)
  signalNames <- grep("RUV3Loess",colnames(dt), value=TRUE)
  #Create a new datatable of annotated normalized signals
  dtn <- dt[,c("Barcode","Well","Ligand","ECMp","Spot","PrintSpot","ArrayRow","ArrayColumn",signalNames), with=FALSE]
  #Add residuals from subtracting the biological medians from each value
  residuals <- dtn[,lapply(.SD,calcResidual), by="Barcode,Well,Ligand,ECMp", .SDcols=signalNames]
  #Add within array location metadata
  residuals$Spot <- as.integer(dtn$Spot)
  residuals$PrintSpot <- as.integer(dt$PrintSpot)
  residuals$ArrayRow <- dtn$ArrayRow
  residuals$ArrayColumn <- dtn$ArrayColumn
  #Create a signal type
  dtn$SignalType <- "Signal"
  residuals$SignalType <- "Residual"
  srDT <- rbind(dtn,residuals)
  
  #Add to carry metadata into matrices
  srDT$BWL <- paste(srDT$Barcode, srDT$Well, srDT$Ligand, sep="_") 
  #srDT$SERC <- paste(srDT$Spot,srDT$ECMp, srDT$ArrayRow, srDT$ArrayColumn, sep="_")
  #srDT$SRC <- paste(srDT$Spot, srDT$ArrayRow, srDT$ArrayColumn, sep="_")
  
  #Set up the M Matrix to denote replicates
  nrControlWells <- sum(grepl("FBS",unique(srDT$BWL[srDT$SignalType=="Signal"])))
  nrLigandWells <- length(unique(srDT$Ligand[srDT$SignalType=="Signal"]))-nrControlWells
  M <-matrix(0, nrow = length(unique(srDT$BWL[srDT$SignalType=="Signal"])), ncol = nrLigandWells+1)
  rownames(M) <- unique(srDT$BWL[srDT$SignalType=="Signal"])
  #Indicate the control wells in the last column
  Mc <- M[grepl("FBS",rownames(M)),]
  Mc[,ncol(Mc)] <-1L
  #Subset to the ligand wells and mark as non-replicate
  Ml <- M[!grepl("FBS",rownames(M)),]
  #Name the rows by their ligand names, removing barcode and well
  MlBWLNames <- rownames(Ml)
  rownames(Ml) <- sub(".*?_","",sub(".*?_","",rownames(Ml)))
  #Name the columns a unique ligand + FBS
  colnames(Ml)<- c(unique(rownames(Ml)),"FBS")
  MlList <- lapply(colnames(Ml),function(cn,Ml){
    #Take the ligand M matrix and put a 1 where the rowname==colname
    Ml[rownames(Ml)==cn,colnames(Ml)==cn] <- 1
    return(Ml[,colnames(Ml)==cn])
  },Ml=Ml)
  Ml <- do.call(cbind, MlList)
  #Restore the row names ot the full BWL values
  rownames(Ml)<- MlBWLNames
  #Add the replicate wells and restore the row order
  M <- rbind(Mc,Ml)
  M <- M[order(rownames(M)),]
  
  srmList <- lapply(signalNames, function(signalName, dt){
    srm <- MEMA:::signalResidualMatrix(dt[,.SD, .SDcols=c("BWL", "PrintSpot", "SignalType", signalName)])
    return(srm)
  },dt=srDT)
  
  names(srmList) <- signalNames
  
  srmRUV3List <- lapply(names(srmList), function(srmName, srmList, M, k){
    Y <- srmList[[srmName]]
    #Hardcode in identification of residuals as the controls
    resStart <- ncol(Y)/2+1
    cIdx=resStart:ncol(Y)
    nY <- RUVIIIArrayWithResiduals(k, Y, M, cIdx, srmName) #Normalize the spot level data
    nY$k <- k
    nY$SignalName <- paste0(srmName,"RUV3")
    setnames(nY,srmName,paste0(srmName,"RUV3"))
    return(nY)
  }, srmList=srmList, M=M, k=k)
  
  #Reannotate with ECMp and MEP
  ECMpDT <- unique(srDT[,list(Well,PrintSpot,Spot,ECMp,ArrayRow,ArrayColumn)])
  
  srmERUV3List <- lapply(srmRUV3List, function(dt,ECMpDT){
    setkey(dt,Well,PrintSpot)
    setkey(ECMpDT,Well,PrintSpot)
    dtECMpDT <- merge(dt,ECMpDT)
    dtECMpDT$MEP <- paste(dtECMpDT$ECMp,dtECMpDT$Ligand,sep="_")
    return(dtECMpDT)
  },ECMpDT=ECMpDT)
  
  
  #Add Loess normalized values for each signal
  RUV3LoessList <- lapply(srmERUV3List, function(dt){
    dtRUV3Loess <- loessNormArray(dt)
  })
  
  #Combine the normalized signal into one data.table
  #with one set of metadata
  signalDT <- do.call(cbind,lapply(RUV3LoessList, function(dt){
    sdt <- dt[,grep("_CP_|_PA_",colnames(dt)), with=FALSE]
  }))
  
  signalMetadataDT <- cbind(RUV3LoessList[[1]][,grep("_CP_|_PA_",colnames(RUV3LoessList[[1]]), invert=TRUE), with=FALSE], signalDT)
  signalMetadataDT <- signalMetadataDT[,SignalName := NULL]
  signalMetadataDT <- signalMetadataDT[,mel := NULL]
  signalMetadataDT <- signalMetadataDT[,Residual := NULL]
  
  #Merge the raw signal back into the dataset
  #signalMetadataDT <- merge(signalMetadataDT,dtn,by=c("Barcode","Well","Spot","ArrayRow","ArrayColumn","ECMp","Ligand"))
  
  return(signalMetadataDT)
}
