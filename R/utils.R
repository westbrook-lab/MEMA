#Read the spot metadata from a non-standard MISVIK gal file
#
#\code{readSpotMetadataMISVIK} returns a dataframe of the contents in a gal
#file and adds array indices.
#
#@param galFile the name of a gal file to be read
#@return A dataframe with the contents of each spot in the gal file. The spots
#  are converted from print block, row and column space to array row and column
#  space. The array space is rotated 180 degrees from the print space. That is,
#  the 1,1,1 position in the print space is the last row and last column of the
#  array space.
#
#  The MISVIK gal file has only 2 print blocks per row.
#
readSpotMetadataMISVIK <- function(galFile){
  #Read the GAL file
  df <- limma::readGAL(galFile)
  layout <- limma::getLayout(df)
  nrCols <- layout$nspot.c*layout$ngrid.c
  nrRows <- layout$nspot.r*layout$ngrid.r
  colnames(df)[colnames(df) == "Block"] <- "Grid"
  #df<-addArrayPositionV2(df,gridsPerRow = layout$ngrid.c)
  df<-addArrayPositionNoRotate(df,gridsPerRow = 2)
  return(df)
}

#Parse names for ECM, Ligand and Concentration values
#
#\code{parseContents} returns a datatable with new columns for the ECM, ligand
#and concentration at each spot.
#
#@param dt a datatable with a Name column to be parsed
#@return The input datatable with new columns for the ECM, Ligand and
#  Concentration values.
#@section Usage: This function assumes a Name column exists and that its
#  format is Concentration_ECM_Ligand or ligand_ECM or Col I_AF 488 or ECM. If Name only
#  contains an ECM, a NULL string is assigned to the Ligand and a
#  concentration value of 0 is assigned to Concentration. If there are 2 words
#  separated by an underscore, the first is the Ligand and the second is the
#  ECM. The concentration is assigned a value of 0. If there are 3 words,
#  the first is assigned to the concentration.
#
#  In the special case of the AF 488 fiducial, the Ligand is assigned the value
#  AF 488.
#
parseContents <- function(dt){
  contentList <- lapply(strsplit(dt$Name,"_"),function(n){
    if(length(n)==1){ECM<-n[[1]]
    Ligand<-""
    Concentration<-0
    return(list(Ligand = Ligand, ECM = ECM,Concentration=Concentration))

    } else if(length(n)==2) {
      Ligand<-n[[1]]
      ECM<-n[[2]]
      if(Ligand=="Col I"&ECM=="AF 488") {
        Ligand<-"AF 488"
        ECM<-"Col I"
      }
      Concentration <- 0
      return(list(Ligand = Ligand, ECM = ECM,Concentration = Concentration))

    } else{
      Concentration <- as.integer(sub("X","",x = n[[1]]))
      Ligand        <- n[[2]]
      ECM           <- n[[3]]
      return(list(Ligand = Ligand, ECM = ECM,Concentration = Concentration))
    }
  }
  )
  contents <- data.table::rbindlist(contentList)
  contents$ConcentrationRank <- match(contents$Concentration,sort(unique(contents$Concentration)))
  dtA <- cbind(dt,contents)
}

#Parse BufferWinners names for ECM, Ligand and Concentration values
#
#\code{parseBufferWinnerContents} returns a datatable with new columns based on the Name column.
#
#@param dt a datatable with a Name column to be parsed
#@return The input datatable with new columns of Condition1, Condition2, Condition3 and Condition4.
#@section Usage: This function assumes a Name column exists and uses an underscore as field separators.
#
parseBufferWinnerContents <- function(dt) {
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    return(list(Condition1 = n[[1]], Condition2 = n[[2]],Condition3 = n[[3]],Condition4 = n[[4]]))
  }
  )
  contents <- data.table::rbindlist(contentList)
  dtA <- cbind(dt,contents)
  return(dtA)
}

#Parse 4wellvalidationContents for ECM, Glycerol, Triton, EDTA and Buffer values
#
#\code{parse4wellvalidationContents} returns a datatable with new columns based on the Name column.
#
#@param dt a datatable with a Name column to be parsed
#@return The input datatable with new columns of ECM, Glycerol, Triton, EDTA and Buffer.
#@section Usage: This functions assumes a Name column exists and uses an underscore as field separators.
#
parse4wellvalidationContents <- function(dt){
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    return(list(ECM = n[[1]], Gycerol = n[[2]],Triton = n[[3]],EDTA = n[[4]],Buffer = n[[5]]))
  }
  )
  contents <- data.table::rbindlist(contentList)
  dtA <- cbind(dt,contents)
  return(dtA)
}

#Parse simulated data for ECM, Ligand and ConcentrationRank values
#
#\code{parseSimulatedContents} returns a datatable with new columns based on
#the Name column.
#
#@param dt a datatable with a Name column to be parsed
#@return The input datatable with new columns of ECM, Ligand and ConcentrationRank.
#@section Usage: This function assumes a Name column exists and uses an
#  underscore as the field separators. The format must be either
#  ECM or ECM_Ligand_ConcentrationRank. ConcentrationRank will be coereced into
#  an integer by removing any non-numeric values.
#
parseSimulatedContents <- function(dt){
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    if (length(n)==2) stop("Invalid format in gal Name field. There should be 1 ECM or 1 ECM with 1 Ligand and a concentration value.")
    if(length(n)==1){ECM<-n[[1]]
    Ligand<-""
    ConcentrationRank<-0
    } else {
      ECM               <- n[[1]]
      Ligand            <- n[[2]]
      ConcentrationRank <- as.integer(gsub("[^[:digit:]]*","",n[[3]]))
    }

    return(list(Ligand = Ligand, ECM = ECM,ConcentrationRank = ConcentrationRank))
  }
  )
  contents <- data.table::rbindlist(contentList)
  dtA <- cbind(dt,contents)
}

#' Coerce the All ECM contents to standard format
#'
#' \code{coerceAllECM}
#'
#'@param df A dataframe read from a gal file
#'@return A dataframe with the Name column reformatted to match the Controlled Vocabulary standard.
#'
#'@export
coerceAllECM <- function(df){
  df$Name <- gsub("__","_",df$Name)
  df$Name <- gsub("PBSCOL I","PBS_none_COL1_UNKNOWN",df$Name)
  df$Name <- gsub("^PBS$","PBS_none",df$Name)
  df$Name <- gsub("COL I$","COL1_UNKNOWN",df$Name)
  df$Name <- gsub("^-$","blank_none",df$Name)

  contentMatrix <- limma::strsplit2(x = df$Name,split="_")

  df$Name <- apply(contentMatrix,1,function(x){
    pUIDs <- paste0(x[seq(1,length(x),2)],"_",x[seq(2,length(x),2)])
    pUIDs <- gsub("^_$","",pUIDs)
    pUID <- paste0(pUIDs,collapse="-")
    pUID <- gsub("[-]*$","",pUID)
  }
  )
  return(df)
}

#' Create a short display version of the LigandAnnotID
#' 
#' @param ligand The name of the ligand
#' @param annotIDs The annotIDs to be simplified
#' @return a vector with the ligand and the last string of each annotID pasted with a _
#' @export
simplifyLigandAnnotID <- function(ligand,annotIDs){
  if(length(annotIDs)){
    splits <- strsplit2(annotIDs, split = "_")
    #Find the last non-empty substring
    us <- apply(splits,1,function(x){
      if(x[length(x)]==""){
        if (length(x)< 2) stop("There are not enough substrings in a ligandAnnotID")
        u <- x[length(x)-1]
      } else{
        u <- x[length(x)]
      }
      return(u)
    })
    ligands <- paste(ligand,us, sep = "_")
  } else ligands <- annotIDs
  return(ligands)
}


#Parse controlled vocabulary data for ECM, Ligand and ConcentrationRank values
#
#\code{parseCVContents} returns a datatable with new columns based on the Name
#column.
#
#@param dt A dataframe or datatable with a Name column to be parsed.
#@return The input datatable with a new column of the ECM names pasted
#  togethers with an underscore.
#@section Usage: This functions assumes a Name column exists and uses a dash
#  symbol - as a field separator between each protein and an underscore between
#  a protein's name and its Uniprot ID.
#
parseCVContents <- function(dt){
  dt$ECM <- lapply(dt$Name,function(x){
    tmp<-strsplit(unlist(strsplit(x,split = "[-]")),split="_")
    ECMNames<-lapply(tmp, function(x){
      x[1]
    })
    paste(ECMNames,collapse="_")
  })
}


# Extract Intensity Endpoints
#
extractIntsEndpts<-function(DT){
  #Use the metadata in the columns Endpoint.xxx to extract the columns of
  #mean intensity data and endpoint name
  epColNames<-grep("(Endpoint)",colnames(DT),value = TRUE)
  intColNames<-paste0("Mean.Intensity.Alexa.",sub("Endpoint.","",epColNames))
  selectColNames<-c(epColNames,intColNames)
  ieDT<-DT[,selectColNames,with=FALSE]
  return(ieDT)
}

# Assign endpoint names
#
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

# Summarize the cell level data to the well level
#
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

# Add the array row, column and index to a GAL file
#
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

