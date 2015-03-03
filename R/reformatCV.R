#Convert the all ECM galfile to a format that uses the pipe separator
#Mark Dane, 3/2/2015

library(limma)
df <- readGAL("~/Documents/MEMATestExp/LI4V008/20150225_LI4V008_AllECMs.gal")

#'Coerce the All ECM contents to standard format
#'
#'\code{coerceALLECM}
#'
#'@param df A dataframe read from a gal file
#'@return A dataframe with the Name column reformatted to match the Controlled Vocabulary standard.
coerceALLECM <- function( df){
  df$Name <- gsub("__","_",df$Name)
  df$Name <- gsub("COL I","COL1_UNKNOWN",df$Name)
  df$Name <- gsub("PBSCOL I","PBS_none_COL1_UNKNOWN",df$Name)
  contentMatrix <- strsplit2(x = df$Name,split="_")

  df$Name <- apply(contentMatrix,1,function(x){
    pUIDs <- paste0(x[seq(1,length(x),2)],"_",x[seq(2,length(x),2)])
    pUIDs <- gsub("^_$","",pUIDs)
    pUID <- paste0(pUIDs,collapse="|")
    pUID <- gsub("[|]*$","",pUID)
  }
  )
  return(df)
}
