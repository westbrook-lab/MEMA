#Normalization functions for processing MEMAs
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
