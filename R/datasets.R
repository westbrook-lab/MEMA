

#' Data from a GAL file and a log file that include spot contents.
#'
#' A dataset containing spotted contents of a MEMA.
#'
#' @format A data frame with 672 rows and 11 variables:
#' \describe{
#'   \item{Row}{print block row}
#'   \item{Column}{print block column}
#'   \item{Block}{print block number}
#'   \item{Name}{name of the spotted contents and if available the uniprot ID separated by an underscore}
#'   \item{ID}{identifier for the source plate well for the spot}
#'   \item{ArrayRow}{array row, starting at the top of the array when well A01 is in the upper left corner}
#'   \item{ArrayColumn}{array column, starting at the left of the array when well A01 is in the upper left corner}
#'   \item{Spot}{an index for the spots that starts at the top left spot and increases from left to right and then from top to bottom}
#'   \item{ShortName}{The part of the Name before the underscore}
#'   \item{PrintOrder}{The order that the spot was printed within the print block}
#'   \item{Depositions}{The number of depositions from the source plate printed at the spot}
#' }
#' @source {Aushon printer}
"spotMetadata"

#' Annotated data from a 4 plate experiment
#'
#' Three channels of annotated population level data from 32 MEMAs.
#' @format A data table with 39 columns and 22,272 rows of data and metadata.
#' @source Tecan LS Reloaded
"popDT"
