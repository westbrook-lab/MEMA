#Functions to clean,merge and normalize data and metadata from one or more well plates.
#Author Mark Dane 11/11/2014
#Updating 12/18/14
#Updating 1/14/2015


#Functions####

#'Read metadata from a multi-sheet excel file
#'
#' \code{readMetadata} reads well level metadata from an excel file. Each sheet in the file represents a different data type. The name of the sheet will become the name of the column. The contents in the sheet will become the data values. Each sheet is foramtted to match a well plate with the A01 well in the upper left corner. First row and column of each sheet are left as labels. The values start in the second row and column.
#' @param xlsFile The name of the excel file.
#' @return a data frame with the data values from each sheet in a column with the name of the sheet.
#' @export
readMetadata<-function(xlsFile){
  #browser()
  sheetList<-sapply(gdata::sheetNames(path.expand(xlsFile)), gdata::read.xls, xls = path.expand(xlsFile), simplify = FALSE,stringsAsFactors=TRUE,check.names=FALSE,row.names="Row/Column")
  nrRows<-dim(sheetList[[1]])[1]
  nrCols<-as.numeric(max(colnames(sheetList[[1]])))
  nrWells=nrRows*nrCols
  sheetDF<-data.frame(lapply(sheetList,function(df,nrCols){
    #create a dataframe from all rows and columns of each sheet
    dfM<-matrix(t(df[,1:nrCols]),byrow=TRUE)
  }, nrCols=nrCols),WellIndex=1:nrWells,Well=wellAN(nrRows,nrCols),check.names=TRUE, stringsAsFactors=FALSE)
  return(sheetDF)
}


#' Merge the Cyto and Main data files
#'
#' \code{mergeCytoMain} merges the Olympus scan^R cytoplasmic and main datatables. The Object.ID in the main datatable is merged with Parent.Object.ID..MO in the cyto datatable.
#' @param cydt The datatable of the cytoplasmic data.
#' @param mdt The datatable of the main data.
#' @return A datatable of the merged data.
#' @export
mergeCytoMain<-function(cydt,mdt){
  #merge the cyto and main data
  #Check that they have the same number of objects
  if(!identical(mdt$Object.ID,cydt$Object.ID)) stop("Main and Cyto files have different Object.ID values and should not be merged.")
  cydt<-data.table::data.table(cydt,key="Parent.Object.ID..MO")
  #Delete unwanted parameters
  keep<-which(!names(cydt) %in% c("Object.ID","Parent.Trace.ID","Well","Parent.Object.ID..Well"))
  cydt<-cydt[,keep,with=FALSE]
  #Get the full names of the intensity and Area columns
  intNames<-grep(pattern="(Intensity|Area)",x=names(cydt),value=TRUE)
  #Add cyto to the parameter names
  data.table::setnames(cydt,intNames,paste0(intNames,"cyto"))
  #Count the number of objects in each well and add as a WellCellCount Column
  mdt<-data.table::data.table(mdt,key="Well")
  wcc<-mdt[,length(Object.ID),by=Well]
  data.table::setnames(wcc,"V1","WellCellCount")
  data.table::setkey(wcc,Well)
  mdt<-mdt[wcc]
  data.table::setkey(mdt,Object.ID)
  #Deleted unwanted parameters
  keep<-which(!names(mdt) %in% c("Parent.Object.ID..MO ","Parent.Trace.ID","Position","Parent.Object.ID..Well"))
  mdt<-mdt[,keep,with=FALSE]
  #merge the cyto and main data files
  merged<-cydt[mdt]
  data.table::setnames(merged,"Parent.Object.ID..MO","Object.ID")
  #If they exist, remove duplicate columns that start with "i."
  select<-grep(pattern = "^(i.)",x = colnames(merged),value=TRUE)
  if(length(select)) merged<-merged[,select:=NULL,with=FALSE]
  return(merged)
}#End mergeCytoMain function

#' Summarize cell data to spot level
#'
#' \code{summarizeToSpot} summarizes cell level data to the spot level.
#' @param cd A datatable of cell level data to be summarized.
#' @return A datatable of the spot level data.
#' @export
summarizeToSpot<-function(cd){
  #browser()
  #Summarize the cell level data to the spot level
  #Create a datatable of the columns to be summariazed
  intNames<-grep(pattern="(Intensity|Area|Spot|WellCellCount|SpotCellCount|^Well$|Objects|EdUPositivePercent|EdUPositiveProportion)",x=names(cd),value=TRUE)
  #Get the metadata column names which are unconstrained
  mdNames<-grep(pattern="(Intensity|Area|WellCellCount|Object.ID|Parent.Object.ID..MO|Elongation.Factor|Circularity.Factor|Perimeter|Max.Feret.Diameter|Center.X|X|Center.Y|Y|Position|SpotCellCount|EdUPositive|EdUCL|EdUPos2N4N)",x=names(cd),value=TRUE,invert=TRUE)
  keep<-cd[,intNames,with=FALSE]
  keepMd<-cd[,mdNames,with=FALSE]
  #Take the mean of each column stratified by spot
  keyCols<-c('Spot','Well')
  wd<-keep[,lapply(.SD,mean),by=keyCols]

  #Get the unique metadata value for each spot
  md<-keepMd[,lapply(.SD,unique),by=keyCols]
  data.table::setkeyv(wd,keyCols)
  data.table::setkeyv(md,keyCols)
  all<-wd[md]
  #Calculate CVs for each set of replicates
  cvNames<-grep(pattern="(Intensity|Area|WellCellCount|SpotCellCount|EdUPositivePercent|EdUPositiveProportion|Name|^Well)",x=names(cd),value=TRUE)
  cvKeep<-all[,cvNames,with=FALSE]
  repSpots<-c('Name','Well')
  cv<-cvKeep[,lapply(.SD,CV),by=repSpots]
  data.table::setnames(cv,colnames(cv), paste0("CV.",colnames(cv)))
  data.table::setkey(cv,CV.Well,CV.Name)
  data.table::setkey(all,Well,Name)
  all<-all[cv]
  return(all)
}


#' Returns a character vector of alphanumeric well names
#'
#' \code{wellAN} is a convenience function that provides alphanumeric names
#' for well plates or arrays.
#' @param nrRows The number of rows in the plate or array.
#' @param nrCols The number of columns in the plate or array.
#' @return A character vector of the well names
#' @export
wellAN<-function(nrRows,nrCols){
  if(nrRows>702)stop("Too many rows to convert. Well alphanumerics are limited to 2 letters")
  Well=unlist(c(lapply(1:nrRows, function(x){
    paste0(paste0(LETTERS[(x-1)%/%26],LETTERS[(1+(x-1)%%26)]), lapply(1:nrCols,function(x){
      if(x<10) x<-paste0("0",x)
      return(as.character(x))
    }))
  })
  ))
  return(Well)
}

#' Change the periods in the column names to spaces
#'
#'\code{cleanColumnNames} substitues a space for any single period in the colmun names of a datatable. This is useful at the end of an anlaysis to have human readable display names.
#' @param dt a datatable
#' @return the same datable with spaces in place of periods.
#' @export
cleanColumnNames<-function(dt){
  data.table::setnames(dt,gsub("[.]"," ",make.names(colnames(dt))))
  return(dt)
}

#' Remove spaces in the column names
#'
#'\code{removeColumnNameSpaces} remove spaces in the colmun names of a datatable. This is useful to compare plates with slight differences in the metadata naming.
#' @param dt a datatable
#' @return the same datable with spaces in the column names removed.
#' @export
removeColumnNameSpaces<-function(dt){
  data.table::setnames(dt,gsub(" ","",make.names(colnames(dt))))
  return(dt)
}

#' merge the data and metadata on the well index column
#'
#' \code{mergeMetadata} merges the well metadata with data on the well value.
#' @param dt A datatable of data values
#' @param mdf A datatable of metadata
#' @return A datatable of annotated data that is merged on the Well values.
#'
#' @export
mergeMetadata<-function(dt, mdf){
  mdt<-data.table::data.table(mdf,key="Well")
  mdt$PrintOrder[is.na(mdt$PrintOrder)] <- 0
  mdt$Depositions[is.na(mdt$Depositions)] <- 0
  #Missing values are read as blanks and not NAs so the check below doesn't help
  if(anyNA(mdt))stop("Some of the metadata is missing.")
  data.table::setkey(dt,"Well")
  merged<-merge(mdt,dt)
}

#' merge the data and metadata on the well index column
#' \code{mergeSpotMetadata} merges the well metadata with data on the well value.
#' @param dt A datatable of data values
#' @param mdf A datatable of metadata
#' @return A datatable of annotated data that is merged on the Well values.
#'
#' @export
mergeSpotMetadata<-function(dt, mdf){
  mdt<-data.table::data.table(mdf,key="Spot")

  #Missing values are read as blanks and not NAs so the check below doesn't help
  if(anyNA(mdt))stop("Some of the metadata is missing.")
  data.table::setkey(dt,"Spot")
  dt$Spot<-(dt$Spot-1)%%(max(mdf$ArrayRow)*max(mdf$ArrayColumn))+1
  merged<-merge(mdt,dt)
}


# remove any wells (and its replicates) that have fewer than wccThresh cells
#
filterWCC<-function(DT,wccThresh){
  #Identify the low cell count wells
  lowWells<-DT$Well[DT$WellCellCount<wccThresh]
  #Identify the replicates of a well
  metaNames<-grep(pattern="(Intensity|Area)",x=names(DT),value=TRUE,invert = TRUE)
  #######start here.....
  #This function uses a WCC column


  return(DT[DT$WellCellCount>wccThresh,])
}

# Normalize on WCC
#
normWCC<-function(DT){
  inChs<-grep("(Alexa.)[[:digit:]]*$|(Alexa.)[[:digit:]]*(cyto10)$",x=names(DT),value = TRUE)
  span=.55
  #Loess normalize the cell data, stratify by channel and cell line
  tmp<-lapply(unique(DT$CellLine),function(cln,DT){
    cDT<-DT[DT$CellLine==cln,]
    normedChs<-lapply(inChs,function(inCh,DT){
      data.table::setkey(DT,WellCellCount)
      set.seed(1234)
      select<-c(1,sample(x = 1:dim(DT)[1],size = dim(DT)[1]*.1,replace = FALSE))
      inChL<-loess(as.formula(paste0(inCh,"~WellCellCount")),data=DT[select,],span=span)
      inChPredicted<-predict(inChL,newdata=DT$WellCellCount)
      return(DT[,inCh,with=FALSE]/inChPredicted)
    },DT=cDT)
    nDT<-do.call("cbind",normedChs)},
    DT=DT)
  #Add the normed versions of the data to the main DT
  data.table::setnames(nDT,names(nDT),paste0(names(nDT),"NormedWCC"))
  nDT<-cbind(DT,nDT)
}


# Add columns of loess normalized SpotCellCount and intensity data
#
addNormalizedInts<-function(DT,span=0.1,...){
  #Get column names of intensity data
  dataNames<-grep(pattern="(Intensity|SpotCellCount)",x=names(DT),value=TRUE)
  #Separate and normalize the data by well
  normedPlateList<-parallel::mclapply(unique(DT$Well), mc.cores=4, function(wellName,DT,span){
    #browser()
    normedWList<-parallel::mclapply(dataNames,mc.cores=4,function(dataName,DTW,span){
      #loess normalize each  column within the current well
      lmodel<-predict(loess(as.formula(paste0(dataName,"~ArrayRow+ArrayColumn")), DTW,span=span))
      normed<-DTW[[dataName]]/lmodel
    },DTW=DT[DT$Well==wellName,],span=span)

    #Combine the normed intensities for each well into a dataframe
    normedDFW<-do.call(data.frame, normedWList)
    colnames(normedDFW)<-paste0('norm.',dataNames)
    #Add the normed intensities to the rest of the data
    DTW<-cbind(DT[DT$Well==wellName,],normedDFW)
  },DT=DT,span=span)
  #
  normedPlate<-do.call(rbind, normedPlateList)
}

# Add columns of median normalized intensity data
#
addMedianNormalizedInts<-function(DT,span=0.1,...){
  #Developed for Teacan population data
  #Get column names of intensity data
  dataNames<-grep(pattern="(Intensity)",x=names(DT),value=TRUE)
  #Separate and normalize the data by well
  normedPlateList<-parallel::mclapply(unique(DT$Well), mc.cores=4, function(wellName,DT,span){
    #browser()
    normedWList<-parallel::mclapply(dataNames,mc.cores=4,function(dataName,DTW,span){
      #loess normalize each  column within the current well
      lmodel<-predict(loess(as.formula(paste0(dataName,"~ArrayRow+ArrayColumn")), DTW,span=span))
      normed<-DTW[[dataName]]/lmodel
    },DTW=DT[DT$Well==wellName,],span=span)

    #Combine the normed intensities for each well into a dataframe
    normedDFW<-do.call(data.frame, normedWList)
    colnames(normedDFW)<-paste0('norm.',dataNames)
    #Add the normed intensities to the rest of the data
    DTW<-cbind(DT[DT$Well==wellName,],normedDFW)
  },DT=DT,span=span)
  #
  normedPlate<-do.call(rbind, normedPlateList)
}

# Extract the named value from a dataset
#
extractNamedValues<-function(normedPlate,namedValue){
  normNames<-grep(pattern="(norm)",x=names(normedPlate),value=TRUE)
  posWList<-lapply(unique(normedPlate$Well),function(well, normedPlate){
    posList<-lapply(normNames,function(normName,DTW){
      posValues<-DTW[DTW$Name==namedValue,normName,with=FALSE]
    },DTW=normedPlate[normedPlate$Well==well,])
    #Bind the positive values for the well into a dataframe
    normedDFW<-do.call(data.frame, posList)
    normedDFW$Well=well
    return(normedDFW)
  },normedPlate=normedPlate)
  #Need to add well metadata
  posW<-do.call(rbind, posWList)
}#End Function

# Add columns of loess normalized intensity data
#
extractNormalizedCntrls<-function(DT,span,...){
  #Get column names of intensity data
  intNames<-grep(pattern="(Intensity)",x=names(DT),value=TRUE)
  #Separate and normalize the data by well
  normedPlateList<-lapply(unique(DT$Well), function(wellName,DT,span){

    normedWList<-lapply(intNames,function(intName,DTW,span){
      #loess normalize each intensity column within the current well
      lmodel<-predict(loess(as.formula(paste0(intName,"~ArrayRow+ArrayColumn")), DTW,span=span))
      normed<-DTW[[intName]]/lmodel
    },DTW=DT[DT$Well==wellName,],span=span)

    #Combine the normed intensities for each well into a dataframe
    normedDFW<-do.call(data.frame, normedWList)
    colnames(normedDFW)<-paste0('norm.',intNames)
    #Add the normed intensities to the rest of the data
    DTW<-cbind(DT[DT$Well==wellName,],normedDFW)
  },DT=DT,span=span)
  #
  normedPlate<-do.call(rbind, normedPlateList)

  #extract the negative and positive normed values for each measurement for each well
  posValues<-data.frame(extractNamedValues(normedPlate,"AllStars Hs Cell Death Control siRNA"),ControlType="Pos")
  negValues<-data.frame(extractNamedValues(normedPlate,"Negative Control siRNA"),ControlType="Neg")
  cntrls<-rbind(posValues, negValues)
}

# Experiment with optimizing on Z prime factor
optimizeZPrimeFactor<-function(DT,...){
  #Get column names of intensity data
  intNames<-grep(pattern="(Intensity)",x=names(DT),value=TRUE)
  #Separate and normalize the data by well
  normedPlateList<-lapply(unique(DT$Well), function(wellName,DT,span){

    normedWList<-lapply(intNames,function(intName,DTW,span){
      #loess normalize each intensity column within the current well
      lmodel<-predict(loess(as.formula(paste0(intName,"~ArrayRow+ArrayColumn")), DTW,span=span))
      normed<-DTW[[intName]]/lmodel
    },DTW=DT[DT$Well==wellName,],span=span)

    #Combine the normed intensities for each well into a dataframe
    normedDFW<-do.call(data.frame, normedWList)
    colnames(normedDFW)<-paste0('norm.',intNames)
    #Add the normed intensities to the rest of the data
    DTW<-cbind(DT[DT$Well==wellName,],normedDFW)
  },DT=DT,span=span)
  #
  normedPlate<-do.call(rbind, normedPlateList)
  tmp<-calculateZPrimeDT(normedPlate)
}

# Calculate the Z prime factor in a datatable
#
calculateZPrimeDT<-function(DT){
  #Given a well of controls values, calculate the ZPrime factors for each channel
  intNames<-grep(pattern="(Intensity)",x=names(DT),value=TRUE)
  ZPrimes<-lapply(intNames, function(intName,DT){
    ZPrime<-calculateZPrime(pos=unlist(DT[DT$ControlType=="Pos",intName, with=FALSE]),
                            neg=unlist(DT[DT$ControlType=="Neg",intName, with=FALSE]))
  },DT=DT)
  ZPrimesDT<-data.frame(ZPrimes)
  colnames(ZPrimesDT)<-intNames
  return(ZPrimesDT)
}

# Calculate Z prime factor
#
calculateZPrime<-function(pos,neg){
  #Given a vector of negative control values and positive control values, calculate
  #the Z Prime Factor value
  1-3*(sd(pos)+sd(neg))/abs(mean(pos)-mean(neg))
}

# Remove unwanted normalized columns
removeNormedColumns<-function(dt){
  #Remove unwanted normalized columns from the dataframe
  normNames<-grep(pattern="(norm)",x=names(dt),value=TRUE)
  if(length(normNames)!=0) dt<-dt[,normNames:=NULL,with=FALSE]
  return(dt)
}


# reorganize the population data by well
#
meltOnWell<-function(DT){
  #browser()
  #DT is a data.table with columns of intensity values separated by
  #barcode, well, wavelength and data type(Raw/Net/Background)
  #This function reorganzies the data into a long format that keeps
  #columns of of intensity readings for each data type and wavelength
  #but moves the barcode and well values into separate columns
  #remove the parts of the name that are specific to this dataset
  intNames<-grep("A1|A2|A3|A4",colnames(DT), value = TRUE)
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

#' Cluster using a Mixture Model
#'
#' \code{autoMMCluster} uses a Mixture Model to assign clusters.
#' @param x A numeric vector
#' @param G the number of groups or clusters
#' @return The cluster assignments of the values in x using Mclust from the mclust package.
#' @export
autoMMCluster<-function(x, G=2)
{# return the cluster assignments for x which is a numeric vector
  #browser()
  x<-data.frame(x)
  cl<-mclust::Mclust(x, G=G,initialization=list(subset=sample(1:nrow(x), size=.05*nrow(x))))[["classification"]]

  return(cl)
}

#' Cluster using kmeans
#'
#' \code{kmeansCluster} is a wrapper function for perfoming kmeans clustering
#' @param x A numeric vector to be clustered
#' @param centers The number of centers or clusters to find.
#' @return The cluster assignments for x using the base kmeans command.
#'
#' @export
kmeansCluster<-function(x, centers=2){
  # return the cluster assignments for x which is a dataframe
  #with the first column as EdU+ clusters and DAPI signal in the 2nd column
  #browser()
  x<-data.frame(x)
  xkmeans<-kmeans(x, centers = centers)
  return(xkmeans[["cluster"]])
}

# Print out correlation plots
#
corPlots <- function(DT,valueName,endPoint){
  browser()
  nrWells <- length(unique(DT$Well))
  if(!nrWells == 4) stop("The number of wells must equal 4")

  #Plot well 1 vs 2
  p <- ggplot(DT, aes_string(x=paste0("DT$",valueName,"[DT$Well== 'A01']"), y=paste0("DT$",valueName,"[DT$Well== 'A02']")))+
    geom_point(size=1)+
    ggtitle(paste(endPoint,"Correlation for PBS in ",unique(DT$Barcode)))
  suppressWarnings(print(p))

  #Plot well 1 vs 3
  p <- ggplot(DT, aes_string(x=paste0("DT$",valueName,"[DT$Well== 'A01']"), y=paste0("DT$",valueName,"[DT$Well== 'A03']")))+
    geom_point(size=1)+
    ggtitle(paste(endPoint,"Correlation for PBS in ",unique(DT$Barcode)))
  suppressWarnings(print(p))

  #Plot well 1 vs 4
  p <- ggplot(DT, aes_string(x=paste0("DT$",valueName,"[DT$Well== 'A01']"), y=paste0("DT$",valueName,"[DT$Well== 'A04']")))+
    geom_point(size=1)+
    ggtitle(paste(endPoint,"Correlation for PBS in ",unique(DT$Barcode)))
  suppressWarnings(print(p))

  #Plot well 2 vs 3
  p <- ggplot(DT, aes_string(x=paste0("DT$",valueName,"[DT$Well== 'A02']"), y=paste0("DT$",valueName,"[DT$Well== 'A03']")))+
    geom_point(size=1)+
    ggtitle(paste(endPoint,"Correlation for PBS in ",unique(DT$Barcode)))
  suppressWarnings(print(p))

  #Plot well 2 vs 4
  p <- ggplot(DT, aes_string(x=paste0("DT$",valueName,"[DT$Well== 'A02']"), y=paste0("DT$",valueName,"[DT$Well== 'A04']")))+
    geom_point(size=1)+
    ggtitle(paste(endPoint,"Correlation for PBS in ",unique(DT$Barcode)))
  suppressWarnings(print(p))

  #Plot well 3 vs 4
  p <- ggplot(DT, aes_string(x=paste0("DT$",valueName,"[DT$Well== 'A03']"), y=paste0("DT$",valueName,"[DT$Well== 'A04']")))+
    geom_point(size=1)+
    ggtitle(paste(endPoint,"Correlation for PBS in ",unique(DT$Barcode)))
  suppressWarnings(print(p))
}


# Deprecated function to support melting population data
castOnBarcodeWell <- function(DT, barcodes){
  wideList <-lapply(barcodes, function(barcode){
    plateList <-lapply(unique(DT$Well[DT$Barcode == barcode]), function(well){
      #subset to only the current barcode and well
      sDT <- DT[DT$Barcode == barcode & DT$Well == well,]
      #Remove the barcode and well
      barcodeWellDT <- sDT[,setnames(.SD, colnames(.SD), paste(barcode,well,colnames(.SD),sep="_")),by="Barcode,Well"]
      #Remove the Barcode and Well Columns
      barcodeWellDT <- barcodeWellDT[,Barcode := NULL]
      barcodeWellDT <- barcodeWellDT[,Well := NULL]
      return(barcodeWellDT)
    })
    #Cbind the DTs in the returned list
    plateDT <- do.call("cbind", plateList)
    return(plateDT)
  })
  #Cbind the DTs in the returned list
  wideDT <- do.call("cbind", wideList)
}

