#Functions to clean,merge and normalize data and metadata from one or more well plates.
#Author Mark Dane 11/11/2014
#Updating 12/18/14
#Updating 1/14/2015


#Functions####

#' Empty function
#'  @export
cleanDT <- function (DT) {
  return(DT)
}

#' Calculate Coefficient of Variation
#'  @export
CV <- function(x){
  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
}

#' Extract Intensity Endpoints
#'  @export
extractIntsEndpts<-function(DT){
  #Use the metadata in the columns Endpoint.xxx to extract the columns of
  #mean intensity data and endpoint name
  epColNames<-grep("(Endpoint)",colnames(DT),value = TRUE)
  intColNames<-paste0("Mean.Intensity.Alexa.",sub("Endpoint.","",epColNames))
  selectColNames<-c(epColNames,intColNames)
  ieDT<-DT[,selectColNames,with=FALSE]
  return(ieDT)
}

#' Assign endpoint names
#'  @export
assignEndPoints<-function(iedt,dt){
  Endpoints<-unique(unlist(iedt[,grep("(Endpoint)",colnames(dt),value = TRUE),with=FALSE]))
  #TODO Generalize this code to determine how many endpoints are in the staining set
  #This code now assumes there are 3 stained endpoints and ignores DAPI
  ei<-do.call("cbind",lapply(Endpoints,function(ep,ieDT){
    #create a column of intensities for each endpoint based on the Endpoint column
    endpointIntsPtrs<-which(ieDT[,1:3,with=FALSE]==ep,arr.ind = TRUE)
    endpointIntsPtrs<-endpointIntsPtrs[order(endpointIntsPtrs[,1]),]
    iValuesM<-as.matrix(ieDT[,4:6,with=FALSE])
    epValues<-iValuesM[endpointIntsPtrs]
  },ieDT=iedt))
  colnames(ei)<-Endpoints
  DT<-cbind(dt,ei)
}

#' Merge the Cyto and Main data files
#' @export
mergeCytoMain<-function(cydt,mdt){
  #merge the cyto and main data
  #browser()
  #Check that they have the same unmber of objects
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

#' Summarize the cell level data to the well level
#' @export
wellLevelData<-function(cd){
  #browser()
  #Create a datatable of the columns to be summariazed
  intNames<-grep(pattern="(Intensity|Area|Well$|WellCellCount)",x=names(cd),value=TRUE)
  #Get the metadata column names which are unconstrained
  if("Spot" %in% names(cd)){
    mdNames<-grep(pattern="(Intensity|Area|WellCellCount|Object.ID|Parent.Object.ID..MO|Elongation.Factor|Circularity.Factor|Perimeter|Max.Feret.Diameter|Center.X|X|Center.Y|Y|Position|WellIndex|Spot|Grid|Column|Row|ID|Name|ArrayRow|ArrayColumn)",x=names(cd),value=TRUE,invert=TRUE)
  } else { mdNames<-grep(pattern="(Intensity|Area|WellCellCount|Object.ID|Parent.Object.ID..MO|Elongation.Factor|Circularity.Factor|Perimeter|Max.Feret.Diameter|Center.X|X|Center.Y|Y|Position|WellIndex)",x=names(cd),value=TRUE,invert=TRUE) }
  keep<-cd[,intNames,with=FALSE]
  keepMd<-cd[,mdNames,with=FALSE]
  #Take the mean of each column stratified by well
  wd<-keep[,lapply(.SD,mean),by=Well]
  #Get the unique metadata value for each well
  #The metadata of a spotted well is not unique
  md<-keepMd[,lapply(.SD,unique),by=Well]
  data.table::setkey(wd,Well)
  data.table::setkey(md,Well)
  all<-wd[md]
  return(all)
}

#' Returns a character vector of alphanumeric well names
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
#' @export
cleanColumnNames<-function(dt){
  data.table::setnames(dt,gsub("[.]"," ",make.names(colnames(dt))))
  return(dt)
}

#' Read metadata from a multi-sheet excel file
#' @export
readMetadata<-function(xlsFile){
  #browser()
  sheetList<-sapply(gdata::sheetNames(xlsFile), gdata::read.xls, xls = xlsFile, simplify = FALSE,stringsAsFactors=TRUE,check.names=FALSE,row.names="Row/Column")
  nrRows<-dim(sheetList[[1]])[1]
  nrCols<-dim(sheetList[[1]])[2]
  nrWells=nrRows*nrCols
  sheetDF<-data.frame(lapply(sheetList,function(df,nrCols){
    #create a dataframe from all rows and columns of each sheet
    dfM<-matrix(t(df[,1:nrCols]),byrow=TRUE)
  }, nrCols=nrCols),WellIndex=1:nrWells,Well=wellAN(nrRows,nrCols),check.names=FALSE, stringsAsFactors=FALSE)
  return(sheetDF)
}

#' merge the data and metadata on the well index column
#' @export
mergeMetadata<-function(dt, mdf){
  mdt<-data.table::data.table(mdf,key="Well")
  #Missing values are read as blanks and not NAs so the check below doesn't help
  if(anyNA(mdt))stop("Some of the metadata is missing.")
  data.table::setkey(dt,"Well")
  merged<-merge(mdt,dt)
}

#' merge the data and metadata on the well index column
#' @export
mergeSpotMetadata<-function(dt, mdf){
  mdt<-data.table::data.table(mdf,key="Spot")
  #Missing values are read as blanks and not NAs so the check below doesn't help
  if(anyNA(mdt))stop("Some of the metadata is missing.")
  data.table::setkey(dt,"Spot")
  dt$Spot<-(dt$Spot-1)%%(max(mdf$ArrayRow)*max(mdf$ArrayColumn))+1
  merged<-merge(mdt,dt)
}

#' Read in a Aushon XML log file to get the deposition counts
#' @export
readLogData<-function(logFile){
  #browser()
  data<-XML::xmlParse(logFile)
  dataList<-XML::xmlToList(data)
  #names(dataList)
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
  return(data)
}

#' remove any wells (and its replicates) that have fewer than wccThresh cells
#' @export
filterWCC<-function(DT,wccThresh){
  #Identify the low cell count wells
  lowWells<-DT$Well[DT$WellCellCount<wccThresh]
  #Identify the replicates of a well
  metaNames<-grep(pattern="(Intensity|Area)",x=names(DT),value=TRUE,invert = TRUE)
  #######start here.....
  #This function uses a WCC column


  return(DT[DT$WellCellCount>wccThresh,])
}

#' Normalize on WCC
#' @export
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

#' Rename the Well column to Spot
#' @export
renameSpottedWells<-function(DT,well){
  data.table::setnames(DT,"Well","Spot")
  DT[,Well:=well]
  return(DT)
}

#' Add columns of loess normalized SpotCellCount and intensity data
#' @export
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

#' Add columns of median normalized intensity data
#' @export
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

#' Extract the named value from a dataset
#' @export
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

#' Add columns of loess normalized intensity data
#' @export
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

#' Experiment with optimizing on Z prime factor
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

#' Calculate the Z prime factor in a datatable
#' @export
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

#' Calculate Z prime factor
#' @export
calculateZPrime<-function(pos,neg){
  #Given a vector of negative control values and positive control values, calculate
  #the Z Prime Factor value
  1-3*(sd(pos)+sd(neg))/abs(mean(pos)-mean(neg))
}

#' Remove unwanted normalized columns
removeNormedColumns<-function(dt){
  #Remove unwanted normalized columns from the dataframe
  normNames<-grep(pattern="(norm)",x=names(dt),value=TRUE)
  if(length(normNames)!=0) dt<-dt[,normNames:=NULL,with=FALSE]
  return(dt)
}

#' Add the array row, column and index to a GAL file
#' @export
addArrayPositionNoRotate<-function(df,gridsPerRow=4){
  #Handle dataframes that have the name block instead of grid
  blockToGrid<-FALSE
  if("Block" %in% colnames(df))
  {
    blockToGrid<-TRUE
    data.table::setnames(df,"Block","Grid")
  }
  #Return the dataframe in the correct space
  #These values are in the printer space
  rowsPerGrid<-max(df$Row)
  colsPerGrid<-max(df$Column)
  #Assign the grid row. These will run from 1:12 when using a 4x12 print head
  df$arrayGridRow<-ceiling(df$Grid/gridsPerRow)
  #Assign the ArrayRow which is in the imaging space
  df$ArrayRow<-(df$arrayGridRow-1)*rowsPerGrid+df$Row
  #Assign the ArrayColumn in imaging space (No rotation)
  df$ArrayColumn<-((df$Grid-1)%%(gridsPerRow))*colsPerGrid+df$Column
  #Rotate the array 180 degrees to match the scan^R
  df$ArrayRow<-max(df$ArrayRow)-df$ArrayRow+1
  df$ArrayColumn<-max(df$ArrayColumn)-df$ArrayColumn+1
  #Order the array by ArrayRow then ArrayColumn
  df<-df[order(df$ArrayRow,df$ArrayColumn),]
  #Remove the arrayGridRow column used for calculations
  df<-subset(df,select=-c(arrayGridRow))
  #Assign a spot number in sequential order by row then column
  df$Spot<-1:(max(df$ArrayColumn)*max(df$ArrayRow))
  #Handle dataframes that have the name block instead of grid
  if(blockToGrid) data.table::setnames(df,"Grid","Block")
  return(df)
}

#' Apply NVS names to staining set
#' @export
nameEndPoints<-function(ssn, dt){
  #rename the intensities into EndPoints and melt the datatable
  if(ssn=="53BP1clPARP"){
    data.table::setnames(dt,"Mean.Intensity.DAPI","DNA")
    data.table::setnames(dt,"Mean.Intensity.Alexa.488","clPARP")
    data.table::setnames(dt,"Mean.Intensity.Alexa.647","p53BP1")
  }else if(ssn=="H3cyclinE"){
    data.table::setnames(dt,"Mean.Intensity.DAPI","DNA")
    data.table::setnames(dt,"Mean.Intensity.Alexa.488","CCNE1")
    data.table::setnames(dt,"Mean.Intensity.Alexa.647","H3meK9")
  } else if(ssn=="Ki67actinITGB1"){
    data.table::setnames(dt,"Mean.Intensity.DAPI","DNA")
    data.table::setnames(dt,"Mean.Intensity.Alexa.488cyto10","Actin")
    data.table::setnames(dt,"Mean.Intensity.Alexa.568","Ki67")
    data.table::setnames(dt,"Mean.Intensity.Alexa.647cyto10","ITGB1")
    #Remove unneeded columns
    rcols<-grep(pattern = "(Mean.Intensity.Alexa.488)|(Mean.Intensity.Alexa.568cyto10)|(Mean.Intensity.Alexa.647)",x = names(dt),value = TRUE)
    dt<-dt[,(rcols):=NULL]
  }else if(ssn=="KRTVIM"){
    data.table::setnames(dt,"Mean.Intensity.DAPI","DNA")
    data.table::setnames(dt,"Mean.Intensity.Alexa.488cyto10","KRT19")
    data.table::setnames(dt,"Mean.Intensity.Alexa.568cyto10","VIM")
    data.table::setnames(dt,"Mean.Intensity.Alexa.647cyto10","KRT14")
    #Remove unneeded columns
    rcols<-grep(pattern = "(Mean.Intensity.Alexa.488)|(Mean.Intensity.Alexa.568)|(Mean.Intensity.Alexa.647)",x = names(dt),value = TRUE)
    dt<-dt[,(rcols):=NULL]
  } else stop("File name has an unknown staining set")
  #Melt the datatable into long format, keeping the id.var columns and putting the
  #EndPoint intensity values into a single column and creating an EndPOint column
  #that contains which Endpoint the intensity is from.
  cdAm<-reshape2::melt(dt,id.vars=c("PlateIndex","drug","biology","ic50_dose_manu","X.ic_20_man","Well","Compound.concentration","Well00" ,"ConcentrationRank", "WellCellCount","Object.ID" ,"Elongation.Factor","Circularity.Factor","Position","Center.Y" ,"X","Y","Center.X","ic50Concentration","Area","Drugset","t_d" ), variable.name="EndPoint", value.name = "Intensity")
  return(cdAm)
}

#' reorganize the population data by well
#' @export
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
#' @export
autoMMCluster<-function(x, G=2)
{# return the cluster assignments for x which is a numeric vector
  #browser()
  x<-data.frame(x)
  cl<-mclust::Mclust(x, G=G,initialization=list(subset=sample(1:nrow(x), size=.05*nrow(x))))[["classification"]]

  return(cl)
}

#' Cluster using kmeans
#' @export
kmeansCluster<-function(x, centers=2)
{# return the cluster assignments for x which is a dataframe
  #with the first column as EdU+ clusters and DAPI signal in the 2nd column
  #browser()
  x<-data.frame(x)
  xkmeans<-kmeans(x, centers = centers)
  return(xkmeans[["cluster"]])
}
