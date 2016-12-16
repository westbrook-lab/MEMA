#Create and test generating the RUV M matrix


#' Create an M matrix for the RUV normalization
#' 
#' The M matrix holds the structure of the dataset in RUV normalization.
#' There is one row for each unit to be normalized and
#' one column for each unique unit type
#' Each row will have one 1 to indicate the unit type
#' all other values will be 0
#' @param dt The datatable to be normalized. There must be a BWL column
#' that defines the unique units, a Ligand column that defines the unique
#' wells and a SignalType column where a value of "Signal" denotes which 
#' BWL and Ligand rows should be included in the M matrix.
#' @return A datatable with values of 1 and 0. 
#' @export
createRUVM <- function(dt)
{
  if(!"BWL" %in% colnames(dt))stop("The data.table to be normalized must have a BWL column")
  if(!"Ligand" %in% colnames(dt))stop("The data.table to be normalized must have a Ligand column")
  if(!"SignalType" %in% colnames(dt))stop("The data.table to be normalized must have a SignalType column")
  #Replace any pipe symbols in the ligand names
  dt$Ligand <- gsub("[/|]","pipe",dt$Ligand)
  dt$BWL <- gsub("[/|]","pipe",dt$BWL)
  #Set up the M Matrix to denote replicate ligand wells
  nrUnits <- length(unique(dt$BWL[dt$SignalType=="Signal"]))
  nrUniqueLigands <- length(unique(dt$Ligand[dt$SignalType=="Signal"]))
  M <-matrix(0, nrow = nrUnits, ncol = nrUniqueLigands)
  rownames(M) <- unique(dt$BWL[dt$SignalType=="Signal"])
  colnames(M) <- unique(dt$Ligand[dt$SignalType=="Signal"])
  gsub(".*_","",(rownames(M)))
  #Indicate the replicate ligands
  for(ligand in colnames(M)){
    #Put a 1 in the rownames that contain the column name
    M[grepl(ligand,rownames(M)),colnames(M)==ligand] <- 1
  }
  #Replace any pipe symbols in the ligand names
  colnames(M) <- gsub("pipe","|",colnames(M))
  rownames(M) <- gsub("pipe","|",rownames(M))
  return(M)
}

#' Create an M matrix for the RUV normalization
#' 
#' The M matrix holds the structure of the dataset in RUV normalization.
#' There is one row for each unit to be normalized and
#' one column for each unique unit type
#' Each row will have one 1 to indicate the unit type
#' all other values will be 0
#' @param dt The datatable to be normalized. There must be a SignalType column 
#' where a value of "Signal" denotes which rows should be included in the M matrix.
#' @param unitID The column name that identifies the names of the units to be normalized.
#' For example, this may have the value BW as the barcode and well can be combined to 
#' create unique identifiers for a unit of data.
#' @param uniqueID The column name that uniquely identifies the replicates. For example,
#' this may have a value of Ligand or LigandDrug depending on the experiment.
#' @return A datatable with values of 1 and 0 that captures the structure of dt.
#' @export
createRUVMGeneral <- function(dt, unitID="BWL", uniqueID="Ligand")
{
  if(!unitID %in% colnames(dt))stop(paste("The data.table to be normalized must have a", unitID, "column"))
  if(!uniqueID %in% colnames(dt))stop(paste("The data.table to be normalized must have a",uniqueID,"column"))
  if(!"SignalType" %in% colnames(dt))stop("The data.table to be normalized must have a SignalType column")
  #Replace any pipe symbols in the ligand names
  dt[[uniqueID]] <- gsub("[/|]","pipe",dt[[uniqueID]])
  dt[[unitID]] <- gsub("[/|]","pipe",dt[[unitID]])
  #Set up the M Matrix to denote replicate ligand wells
  nrUnits <- length(unique(dt[[unitID]][dt$SignalType=="Signal"]))
  nrUniqueIDs <- length(unique(dt[[uniqueID]][dt$SignalType=="Signal"]))
  M <-matrix(0, nrow = nrUnits, ncol = nrUniqueIDs)
  rownames(M) <- unique(dt[[unitID]][dt$SignalType=="Signal"])
  colnames(M) <- unique(dt[[uniqueID]][dt$SignalType=="Signal"])
  #Indicate the replicate ligands
  for(uID in rownames(M)){
    #For each row, put a 1 in column that matches it's uniqueID value
    M[uID,dt[[uniqueID]][dt[[unitID]]==uID]==colnames(M)] <- 1
  }
  #Replace any pipe symbols in the ligand names
  colnames(M) <- gsub("pipe","|",colnames(M))
  rownames(M) <- gsub("pipe","|",rownames(M))
  return(M)
}

