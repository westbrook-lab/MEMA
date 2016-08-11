#Normalization functions for processing MEMAs

#' Median normalize a vector
#' 
#' @param DT A datatable with a column \code{value}
#' @param value The name of the column in DT to be normalized
#' 
#' @return A numeric vector that is \code{DT[[value]]} dvided by its median
#' @export
#' 
medianNorm <- function(DT, value){
  normedValues <- DT[, value, with = FALSE]/median(unlist(DT[, value, with = FALSE]), na.rm=TRUE)
}

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
#'   normalized, a \code{ECMp} column that has the printed ECM names and a
#'   \code{Growth.Factors} or \code{Ligand}column that has the growth factor names.
#' @param value The name of the column of values to be normalized
#' @param baseECM A regular expression for the name or names of the printed ECM(s) to be normalized against
#' @param baseGF A regular expression for the name or names of the soluble growth factors to be normalized against
#' @return A numeric vector of the normalized values
#'
#' @section Details: \code{normWellsWithinPlate} normalizes the value column of
#'   all MEPs by dividing the median value of the replicates of the MEP that
#'   is the pairing of baseECM  with baseGF.
#'   @export
normWellsWithinPlate <- function(DT, value, baseECM, baseGF) {
  if(!c("ECMp") %in% colnames(DT)) stop(paste("DT must contain a ECMp column."))
  if(!c(value) %in% colnames(DT)) stop(paste("DT must contain a", value, "column."))
  if("Ligand" %in% colnames(DT)){
    valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ECMp)  & grepl(baseGF,DT$Ligand)),value, with=FALSE]), na.rm = TRUE)
  } else if (c("Growth.Factors") %in% colnames(DT)) {
    valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ECMp)  & grepl(baseGF,DT$Growth.Factors)),value, with=FALSE]), na.rm = TRUE)
  } else stop (paste("DT must contain a Growth.Factors or Ligand column."))
  normedValues <- DT[,value,with=FALSE]/valueMedian
  return(normedValues)
}

#' Create a median normalized loess model of an array
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

#' RZS Normalize a Column of Data
#'
#' \code{normRZSWellsWithinPlate} normalizes all elements of DT[[value]] by
#' subtracting the median of DT[[value]] of all baseECM spots in the
#' baseL wells, then divides the result by the MAD*1.48 of all baseECM spots in
#' the baseL wells
#'@param DT A datatable with value, baseECM and baseL, ECMp and
#'Ligand columns
#'@param value A single column name of the value to be normalized
#'@param baseECM A single character string or a regular expression that selects
#'the ECM(s) that are used as the base for normalization.
#'@param baseL A single character string or a regular expression that selects
#'the ligand used as the base for normalization.
#'@return a vector of RZS normalized values
#' @export
#'
normRZSWellsWithinPlate <- function(DT, value, baseECM, baseL) {
  if(!"ECMp" %in% colnames(DT)) stop (paste("DT must contain an ECMp column."))
  if(!"Ligand" %in% colnames(DT)) stop (paste("DT must contain a Ligand column."))
  if(!c(value) %in% colnames(DT)) stop(paste("DT must contain a", value, "column."))

  valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ECMp) & grepl(baseL,DT$Ligand)), value, with=FALSE]), na.rm = TRUE)
  if (is.na(valueMedian)) stop(paste("Normalization calculated an NA median for",value, baseECM, baseL))

  valueMAD <- mad(unlist(DT[(grepl(baseECM, DT$ECMp)  & grepl(baseL,DT$Ligand)),value, with=FALSE]), na.rm = TRUE)
  #Correct for 0 MAD values
  valueMAD <- valueMAD+.01
  normedValues <- (DT[,value,with=FALSE]-valueMedian)/valueMAD
  return(normedValues)
}

#' Normalize selected values in a dataset on a plate basis
#' 
#' A wrapper function for \code{normRZSWellsWithinPlate} that selects the
#' \code{_CP_|_QI_|_PA_|SpotCellCount|Lineage} columns of dt if they exist and 
#' normalizes them on a plate basis
#' @param dt A data.table with a \code{Barcode} column numeric values to be RZS normalized 
#'  using all ECM proteins in the FBS well
#' @return A datatable with the normalized values
#' @export
normRZSDataset <- function(dt){
  parmNormedList <- lapply(grep("_CP_|_QI_|_PA_|SpotCellCount|Lineage",colnames(dt),value = TRUE), function(parm){
    dt <- dt[,paste0(parm,"_RZSNorm") := normRZSWellsWithinPlate(.SD, value=parm, baseECM = ".*",baseL = "FBS"), by="Barcode"]
    return(dt)
  })
  return(parmNormedList[[length(parmNormedList)]])
}

#'Apply RUV3 and Loess Normalization to the signals in a dataset
#' @export
normRUV3LoessResidualsDisplay <- function(dt, k){
  setkey(dt,Barcode,Well,Ligand,ECMp)
  metadataNames <- "Barcode|Well|^Spot$|^PrintSpot$|ArrayRow|ArrayColumn|^ECMp$|^Ligand$"
  signalNames <- grep(metadataNames,colnames(dt),invert=TRUE, value=TRUE)
  
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
  #srDT$WsRC <- paste(srDT$PrintSpot, srDT$ArrayRow, srDT$ArrayColumn, sep="_")
  
  
  #Set up the M Matrix to denote replicates
  nrControlWells <- sum(grepl("FBS",unique(srDT$BWL[srDT$SignalType=="Signal"])))
  nrLigandWells <- length(unique(srDT$BWL[srDT$SignalType=="Signal"]))-nrControlWells
  M <-matrix(0, nrow = length(unique(srDT$BWL[srDT$SignalType=="Signal"])), ncol = nrLigandWells+1)
  rownames(M) <- unique(srDT$BWL[srDT$SignalType=="Signal"])
  #Indicate the control wells in the last column
  Mc <- M[grepl("FBS",rownames(M)),]
  Mc[,ncol(Mc)] <-1L
  #Subset to the ligand wells and mark as non-replicate
  Ml <- M[!grepl("FBS",rownames(M)),]
  for(i in 1:nrLigandWells) {
    Ml[i,i] <- 1
  }
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
    nY <- RUVIIIArrayWithResiduals(k, Y, M, cIdx, srmName, verboseDisplay = TRUE)[["nYm"]] #Normalize the spot level data
    nY$k <- k
    nY$SignalName <- paste0(srmName,"RUV3")
    setnames(nY,srmName,paste0(srmName,"RUV3"))
    nY[[srmName]] <- as.vector(Y[,1:(resStart-1)])
    return(nY)
  }, srmList=srmList, M=M, k=k)
  
  #Reannotate with ECMp, MEP, ArrayRow and ArrayColumn
  ECMpDT <- unique(srDT[,list(Well,PrintSpot,Spot,ECMp,ArrayRow,ArrayColumn)])
  
  srmERUV3List <- lapply(srmRUV3List, function(dt,ECMpDT){
    setkey(dt,Well,PrintSpot)
    setkey(ECMpDT,Well,PrintSpot)
    dtECMpDT <- merge(dt,ECMpDT)
    dtECMpDT$MEP <- paste(dtECMpDT$ECMp,dtECMpDT$Ligand,sep="_")
    return(dtECMpDT)
  },ECMpDT=ECMpDT)
  
  
  #Add Loess normalized values for each RUV3 normalized signal
  RUV3LoessList <- lapply(srmERUV3List, function(dt){
    dtRUV3Loess <- loessNormArray(dt)
  })
  
  #Add Loess normalized values for each Raw signal
  RUV3LoessList <- lapply(srmERUV3List, function(dt){
    dt$SignalName <- sub("RUV3","",dt$SignalName)
    dtLoess <- loessNormArray(dt)
  })
  
  #Combine the normalized signal into one data.table
  #with one set of metadata
  signalDT <- do.call(cbind,lapply(RUV3LoessList, function(dt){
    sdt <- dt[,grep("_CP_|_PA_|Cells|Reference",colnames(dt)), with=FALSE]
  }))
  
  signalMetadataDT <- cbind(RUV3LoessList[[1]][,grep("_CP_|_PA_|Cells|Reference",colnames(RUV3LoessList[[1]]), invert=TRUE), with=FALSE], signalDT)
  signalMetadataDT <- signalMetadataDT[,SignalName := NULL]
  signalMetadataDT <- signalMetadataDT[,mel := NULL]
  signalMetadataDT <- signalMetadataDT[,Residual := NULL]
  return(signalMetadataDT)
}

#'Apply RUV3 and Loess Normalization to the signals in a dataset
#' @export
normRUV3LoessResiduals <- function(dt, k){
  setkey(dt,Barcode,Well,Ligand,ECMp)
  metadataNames <- "Barcode|Well|^Spot$|^PrintSpot$|ArrayRow|ArrayColumn|^ECMp$|^Ligand$"
  signalNames <- grep(metadataNames,colnames(dt),invert=TRUE, value=TRUE)
  
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
  #srDT$WsRC <- paste(srDT$PrintSpot, srDT$ArrayRow, srDT$ArrayColumn, sep="_")
  
  
  #Set up the M Matrix to denote replicates
  nrControlWells <- sum(grepl("FBS",unique(srDT$BWL[srDT$SignalType=="Signal"])))
  nrLigandWells <- length(unique(srDT$BWL[srDT$SignalType=="Signal"]))-nrControlWells
  M <-matrix(0, nrow = length(unique(srDT$BWL[srDT$SignalType=="Signal"])), ncol = nrLigandWells+1)
  rownames(M) <- unique(srDT$BWL[srDT$SignalType=="Signal"])
  #Indicate the control wells in the last column
  Mc <- M[grepl("FBS",rownames(M)),]
  Mc[,ncol(Mc)] <-1L
  #Subset to the ligand wells and mark as non-replicate
  Ml <- M[!grepl("FBS",rownames(M)),]
  for(i in 1:nrLigandWells) {
    Ml[i,i] <- 1
  }
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
    nY[[srmName]] <- as.vector(Y[,1:(resStart-1)])
    return(nY)
  }, srmList=srmList, M=M, k=k)
  
  #Reannotate with ECMp, MEP, ArrayRow and ArrayColumn
  ECMpDT <- unique(srDT[,list(Well,PrintSpot,Spot,ECMp,ArrayRow,ArrayColumn)])
  
  srmERUV3List <- lapply(srmRUV3List, function(dt,ECMpDT){
    setkey(dt,Well,PrintSpot)
    setkey(ECMpDT,Well,PrintSpot)
    dtECMpDT <- merge(dt,ECMpDT)
    dtECMpDT$MEP <- paste(dtECMpDT$ECMp,dtECMpDT$Ligand,sep="_")
    return(dtECMpDT)
  },ECMpDT=ECMpDT)
  
  
  #Add Loess normalized values for each RUV3 normalized signal
  RUV3LoessList <- lapply(srmERUV3List, function(dt){
    dtRUV3Loess <- loessNormArray(dt)
  })
  
  #Add Loess normalized values for each Raw signal
  RUV3LoessList <- lapply(srmERUV3List, function(dt){
    dt$SignalName <- sub("RUV3","",dt$SignalName)
    dtLoess <- loessNormArray(dt)
  })
  
  #Combine the normalized signal into one data.table
  #with one set of metadata
  signalDT <- do.call(cbind,lapply(RUV3LoessList, function(dt){
    sdt <- dt[,grep("_CP_|_PA_|Cells|Reference",colnames(dt)), with=FALSE]
  }))
  
  signalMetadataDT <- cbind(RUV3LoessList[[1]][,grep("_CP_|_PA_|Cells|Reference",colnames(RUV3LoessList[[1]]), invert=TRUE), with=FALSE], signalDT)
  signalMetadataDT <- signalMetadataDT[,SignalName := NULL]
  signalMetadataDT <- signalMetadataDT[,mel := NULL]
  signalMetadataDT <- signalMetadataDT[,Residual := NULL]
  return(signalMetadataDT)
}


#' Loess normalize an array using the spatial residuals
loessNorm <- function(Value,Residual,ArrayRow,ArrayColumn){
  dt <-data.table(Value=Value,Residual=Residual,ArrayRow=ArrayRow,ArrayColumn=ArrayColumn)
  lm <- loess(Residual~ArrayRow+ArrayColumn, dt, span=.7)
  dt$ResidualLoess<-predict(lm)
  dt <- dt[,ValueLoess := Value-ResidualLoess]
  return(ValueLoess = dt$ValueLoess)
}

#' Loess normalize values within an array
#'@export
loessNormArray <- function(dt){
  #Identify the Signal name
  signalName <- unique(dt$SignalName)
  setnames(dt,signalName,"Value")
  #Get the median of the replicates within the array
  dt <- dt[,mel := median(Value), by=c("BWL","ECMp")]
  #Get the residuals from the spot median
  dt <- dt[,Residual := Value-mel]
  #Subtract the loess model of each array's residuals from the signal
  dt <- dt[, ValueLoess:= loessNorm(Value,Residual,ArrayRow,ArrayColumn), by="BWL"]
  setnames(dt,"Value", signalName)
  setnames(dt,"ValueLoess", paste0(signalName,"Loess"))
  return(dt)
}


#' Calculate the residuals from the median of a vector of numeric values
#' @export
calcResidual <- function(x){
  mel <- as.numeric(median(x, na.rm=TRUE))
  return(x-mel)
}

#' Apply RUV3 normalization on a signal and its residuals
#' 
#' Assumes there are signal values in the first half of each row
#' and residuals in the second half
#' @export
RUVIIIArrayWithResiduals <- function(k, Y, M, cIdx, signalName, verboseDisplay=FALSE){
  YRUVIII <- RUVIII(Y, M, cIdx, k)
  nY <- YRUVIII[["newY"]]
  #Remove residuals
  nY <- nY[,1:(ncol(nY)/2)]
  #melt matrix to have Spot and Ligand columns
  nYm <- melt(nY, varnames=c("BWL","PrintSpot"), value.name=signalName)
  nYm$BWL <- as.character(nYm$BWL)
  #nYm$SRC <- as.character(nYm$SRC)
  splits <- limma::strsplit2(nYm$BWL,split = "_")
  nYm$Barcode <- splits[,1]
  nYm$Well <- splits[,2]
  nYm$Ligand <- sub(".*?_","",sub(".*?_","",nYm$BWL))
  if(verboseDisplay){
    return(list(nYm=data.table(nYm), fullAlpha=YRUVIII[["fullalpha"]], W=RUVIII[["W"]]))
  }
  return(data.table(nYm))
}

#' Generate a matrix of signals and residuals
#' 
#' This function reorganizes a data.table into a matrix suitable for
#' RUV3 normalization.
#' 
#'@param dt a data.table with columns named BWL, SRC, SignalType and a column 
#' named by the unique value in the SignalType column.
#'@return A numeric matrix with BWL rows and two sets of SERC columns. The second
#'set of SRC columns are the residuals from the medians of each column and have
#'"_Rel" appended to their names.
signalResidualMatrix <- function(dt){
  signalName <- colnames(dt)[ncol(dt)]
  if(grepl("Logit", signalName)){
    fill <- log2(.01/(1-.01))
  } else if(grepl("Log", signalName)){
    fill <- log2(.001)
  } else {
    fill <- 0
  }

  dts <- data.table(dcast(dt[dt$SignalType=="Signal",], BWL~PrintSpot, value.var=signalName, fill=fill, na.rm=TRUE))
  dtr <- data.table(dcast(dt[dt$SignalType=="Residual",], BWL~PrintSpot, value.var=signalName, fill=fill, na.rm=TRUE))
  rowNames <- dts$BWL
  dts <- dts[,BWL := NULL]
  dtr <- dtr[,BWL:=NULL]
  setnames(dtr,colnames(dtr),paste0(colnames(dtr),"_Rel"))
  dtsr <- cbind(dts,dtr)
  srm <- matrix(unlist(dtsr),nrow=nrow(dtsr))
  rownames(srm) <- rowNames
  colnames(srm) <- colnames(dtsr)
  return(srm)
}

#' Perform RUV3 removal of unwanted variations
#' This function is written by Johann Gagnon-Bartsch and will become part of the
#' ruv package.
RUVIII = function(Y, M, ctl, k=NULL, eta=NULL, average=FALSE, fullalpha=NULL)
{
  Y = RUV1(Y,eta,ctl)
  if (is.null(k))
  {
    ycyctinv = solve(Y[,ctl]%*%t(Y[,ctl]))
    newY = (M%*%solve(t(M)%*%ycyctinv%*%M)%*%(t(M)%*%ycyctinv)) %*% Y
    fullalpha=NULL
  }
  else if (k == 0) newY = Y
  else
  {
    m = nrow(Y)
    Y0 = residop(Y,M)
    fullalpha = t(svd(Y0%*%t(Y0))$u[,1:(m-ncol(M)),drop=FALSE])%*%Y
    k<-min(k,nrow(fullalpha))
    alpha = fullalpha[1:k,,drop=FALSE]
    ac = alpha[,ctl,drop=FALSE]
    W = Y[,ctl]%*%t(ac)%*%solve(ac%*%t(ac))
    newY = Y - W%*%alpha
  }
  if (average) newY = ((1/apply(M,2,sum))*t(M)) %*% newY
  return(list(newY = newY, fullalpha=fullalpha, W=W))
}

#' Apply the RUV3 algortihm 
normRUV3Dataset <- function(dt, k){
  #browser()
  #Setup data with plate as the unit
  #There are 694 negative controls and all plates are replicates
  
  setkey(dt, Barcode,Well,Spot)     #Sort the data
  metadataNames <- "Barcode|Well|^Spot$"
  signalNames <- grep(metadataNames,colnames(dt),invert=TRUE, value=TRUE)
  
  dt$WS <- paste(dt$Well, dt$Spot,sep="_") #Add to carry metadata to matrix
  
  nYL <- lapply(signalNames, function(signal, dt, M){
    #Create appropriate fill for missing values for each signal
    if(grepl("EdU|Proportion|Ecc", signal)){
      fill <- log2(.01/(1-.01))
    } else if(grepl("Log", signal)){
      fill <- log2(.001)
    } else {
      fill <- 0
    }
    
    #Cast into barcode rows and well spot columns
    dtc <- dcast(dt, Barcode~WS, value.var=signal, fill=fill)
    #Remove the Barcode column and use it as rownames in the matrix
    barcodes <-dtc$Barcode
    dtc <- dtc[,Barcode := NULL]
    Y <- matrix(unlist(dtc), nrow=nrow(dtc), dimnames=list(barcodes, colnames(dtc)))
    k<-min(k, nrow(Y)-1)
    cIdx <- which(grepl("A03",colnames(Y)))
    nY <- MEMA:::RUVIII(Y, M, cIdx, k)[["newY"]]
    #melt matrix to have ECMp and Ligand columns
    nYm <- melt(nY, varnames=c("Barcode","WS"),  as.is=TRUE)
    nYm <- data.table(nYm)
    
    #Add the name of the signal and convert back to well and spot
    nYm$Signal <- paste0(signal,"_Norm")
    return(nYm)
    
  }, dt=dt, M=matrix(1,nrow=length(unique(dt$Barcode)))
  # ,mc.cores = detectCores())
  )
  
  nYdtmelt <- rbindlist(nYL)
  nY <- dcast(nYdtmelt, Barcode+WS~Signal, value.var="value")
  nY$Well <- gsub("_.*","",nY$WS)
  nY$Spot <- as.integer(gsub(".*_","",nY$WS))
  nY <- nY[,WS:=NULL]
  return(nY)
}

