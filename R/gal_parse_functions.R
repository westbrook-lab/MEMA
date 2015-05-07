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
#'  space. The array space is aligned with the print space. That is,
#'  the 1,1,1 position in the print space is the first row and first column of the
#'  array space.
#' @export
readSpotMetadata <- function(galFile) {
  #Read the GAL file
  #browser()
  df <- limma::readGAL(galFile)
  layout <- limma::getLayout(df)
  nrCols <- layout$nspot.c*layout$ngrid.c
  nrRows <- layout$nspot.r*layout$ngrid.r
  colnames(df)[colnames(df) == "Grid"] <- "Block"
  df<-addArrayPositionNoRotate(df,gridsPerRow = layout$ngrid.c)
  DT<-data.table::data.table(df,key=c("Block","Row","Column"))
  return(DT)
}

