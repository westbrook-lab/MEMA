#Functions to parse gal files
#Author Mark Dane 3/2/2015

#'Read the spot metadata from a gal file
#'
#'\code{readSpotMetadata} returns a dataframe of the contents in a gal file and
#'adds array indices.
#'
#'@param galFile the name of a gal file to be read
#'@return A dataframe with the contents of each spot in the gal file. The spots
#'  are converted from print block, row and column space to array row and column
#'  space. The array space is rotated 180 degrees from the print space. That is,
#'  the 1,1,1 position in the print space is the last row and last column of the
#'  array space.
readSpotMetadata <- function(galFile) {
  #Read the GAL file
  #browser()
  df <- readGAL(galFile)
  layout <- getLayout(df)
  nrCols <- layout$nspot.c*layout$ngrid.c
  nrRows <- layout$nspot.r*layout$ngrid.r
  colnames(df)[colnames(df) == "Block"] <- "Grid"
  df<-addArrayPositionNoRotate(df,gridsPerRow = layout$ngrid.c)
  return(df)
}

#'Read the spot metadata from a non-standard MISVIK gal file
#'
#'\code{readSpotMetadataMISVIK} returns a dataframe of the contents in a gal
#'file and adds array indices.
#'
#'@param galFile the name of a gal file to be read
#'@return A dataframe with the contents of each spot in the gal file. The spots
#'  are converted from print block, row and column space to array row and column
#'  space. The array space is rotated 180 degrees from the print space. That is,
#'  the 1,1,1 position in the print space is the last row and last column of the
#'  array space.
#'
#'  The MISVIK gal file has only 2 print blocks per row.
readSpotMetadataMISVIK <- function(galFile){
  #Read the GAL file
  df <- readGAL(galFile)
  layout <- getLayout(df)
  nrCols <- layout$nspot.c*layout$ngrid.c
  nrRows <- layout$nspot.r*layout$ngrid.r
  colnames(df)[colnames(df) == "Block"] <- "Grid"
  #df<-addArrayPositionV2(df,gridsPerRow = layout$ngrid.c)
  df<-addArrayPositionNoRotate(df,gridsPerRow = 2)
  return(df)
}

#'Parse names for ECM, Ligand and Concentration values
#'
#'\code{parseContents} returns a datatable with new columns for the ECM, ligand
#'and concentration at each spot.
#'
#'@param dt a datatable with a Name column to be parsed
#'@return The input datatable with new columns for the ECM, Ligand and
#'  Concentration values.
#'@section Usage: This functions assumes a Name column exists and that its
#'  format is ECM_Ligand_Concentration. If there are no underscores, then the
#'  Name is interpreted as an ECM, a NULL string is assigned to the Ligand and a
#'  concentration value of 0 is assigned to Concentration. If there are 2 words
#'  separated by an underscore, the first is the ECM and the second is the
#'  Ligand. The concentration is assigned a value of 0. If there are 3 words,
#'  the first is assigned to the concentration.
#'
#'  In the special case of the AF 488 fiducial, the Ligand is assigned the value
#'  AF 488.
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
  contents <- rbindlist(contentList)
  contents$ConcentrationRank <- match(contents$Concentration,sort(unique(contents$Concentration)))
  dtA <- cbind(dt,contents)
}

#'Parse BufferWinners names for ECM, Ligand and Concentration values
#'
#'\code{parseBufferWinnerContents} returns a datatable with new columns based on the Name column.
#'
#'@param dt a datatable with a Name column to be parsed
#'@return The input datatable with new columns of Condition1, Condition2, Condition3 and Condition4.
#'@section Usage: This functions assumes a Name column exists uses an underscore as field separators.
parseBufferWinnerContents <- function(dt) {
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    return(list(Condition1 = n[[1]], Condition2 = n[[2]],Condition3 = n[[3]],Condition4 = n[[4]]))
  }
  )
  contents <- rbindlist(contentList)
  dtA <- cbind(dt,contents)
  return(dtA)
}

#'Parse 4wellvalidationContents for ECM, Glycerol, Triton, EDTA and Buffer values
#'
#'\code{parse4wellvalidationContents} returns a datatable with new columns based on the Name column.
#'
#'@param dt a datatable with a Name column to be parsed
#'@return The input datatable with new columns of ECM, Glycerol, Triton, EDTA and Buffer.
#'@section Usage: This functions assumes a Name column exists and uses an underscore as field separators.
parse4wellvalidationContents <- function(dt){
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    return(list(ECM = n[[1]], Gycerol = n[[2]],Triton = n[[3]],EDTA = n[[4]],Buffer = n[[5]]))
  }
  )
  contents <- rbindlist(contentList)
  dtA <- cbind(dt,contents)
  return(dtA)
}

#'Parse simulated data for ECM, Ligand and ConcentrationRank values
#'
#'\code{parseSimulatedContents} returns a datatable with new columns based on
#'the Name column.
#'
#'@param dt a datatable with a Name column to be parsed
#'@return The input datatable with new columns of ECM, Ligand and ConcentrationRank.
#'@section Usage: This functions assumes a Name column exists and uses an
#'  underscore as field separators. The format is assumed to be
#'  ECM_Ligand_ConcentrationRank
parseSimulatedContents <- function(dt){
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    ECM               <- n[[1]]
    Ligand            <- n[[2]]
    ConcentrationRank <-as.integer(n[[3]])
    return(list(Ligand = Ligand, ECM = ECM,ConcentrationRank = ConcentrationRank))
  }

  )
  contents <- rbindlist(contentList)
  dtA <- cbind(dt,contents)
}

#'Coerce the All ECM contents to standard format
#'
#'\code{coerceAllECM}
#'
#'@param df A dataframe read from a gal file
#'@return A dataframe with the Name column reformatted to match the Controlled Vocabulary standard.
coerceAllECM <- function(df){
  df$Name <- gsub("__","_",df$Name)
  df$Name <- gsub("PBSCOL I","PBS_none_COL1_UNKNOWN",df$Name)
  df$Name <- gsub("^PBS$","PBS_none",df$Name)
  df$Name <- gsub("COL I$","COL1_UNKNOWN",df$Name)
  df$Name <- gsub("^-$","blank_none",df$Name)

  contentMatrix <- strsplit2(x = df$Name,split="_")

  df$Name <- apply(contentMatrix,1,function(x){
    pUIDs <- paste0(x[seq(1,length(x),2)],"_",x[seq(2,length(x),2)])
    pUIDs <- gsub("^_$","",pUIDs)
    pUID <- paste0(pUIDs,collapse="-")
    pUID <- gsub("[-]*$","",pUID)
  }
  )
  return(df)
  }

#'Parse controlled vocabulary data for ECM, Ligand and ConcentrationRank values
#'
#'\code{parseCVContents} returns a datatable with new columns based on the Name
#'column.
#'
#'@param dt A dataframe or datatable with a Name column to be parsed
#'@return The input datatable with a new column of the ECM names pasted
#'  togethers with an underscore.
#'@section Usage: This functions assumes a Name column exists and uses a pipe
#'  symbol | as a field separator between each protein and an underscore between
#'  a protein's name and its Uniprot ID.
parseCVContents <- function(dt){
  dt$ECM <- lapply(dt$Name,function(x){
    tmp<-strsplit(unlist(strsplit(x,split = "[|]")),split="_")
    ECMNames<-lapply(tmp, function(x){
      x[1]
    })
    paste(ECMNames,collapse="_")
  })
}

