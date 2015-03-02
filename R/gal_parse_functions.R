#Functions to parse gal files
#Author Mark Dane 3/2/2015

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

parseContents <- function(dt){
  #parse the content names in the gal file for Ligand, ECM and concentrations
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

parse4wellvalidationContents <- function(dt){
  #browser()
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    return(list(ECM = n[[1]], Gycerol = n[[2]],Triton = n[[3]],EDTA = n[[4]],Buffer = n[[5]]))
  }
  )
  contents <- rbindlist(contentList)
  dtA <- cbind(dt,contents)
  return(dtA)
}

parseSimulatedContents <- function(dt){
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    ECM               <- n[[1]]
    Ligand            <-n [[2]]
    ConcentrationRank <-as.integer(n[[3]])
    return(list(Ligand = Ligand, ECM = ECM,ConcentrationRank = ConcentrationRank))
  }

  )
  contents <- rbindlist(contentList)
  dtA <- cbind(dt,contents)
}

parseCVContents <- function(dt){
  #parse a gal file that uses the controlled vocabulary format
  #Valid formats are:
  #ECMName_UniprotID
  #ECMName_UniprotID_ECMName_UniprotID an inlimited number of times
  #ECMName_UniprotID|LigandName_UniprotID
  #an unlimited number of ECMs and Ligands

  contentList <- lapply(strsplit(dt$Name,"|"))


}
