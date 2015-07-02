#Functions to support MEMAs printed in 8 well plates

#'\code{spotCellDensities} Calculate the neghborhood density around each cell
#'@param spot A datatable with X and Y columns
#'@param radius The radial distance that defines the neighborhood around the cell
#'@return A numeric vector of length nrow(spot) with the cell density values
spotCellDensities <- function(spot,radius=(max(spot$X)-min(spot$X))/5) {
  distMatrix <- as.matrix(dist(spot[,list(X,Y)]))
  count <- apply(distMatrix, 2, function(x){sum(x <= radius) - 1})
  cellDensity <- count/(pi*radius^2)
  return(cellDensity)
}


#' Calculate positional derived parameters for cell-level data
#'@param DT A cell level datatable with plate-level x,y coordinates
#'@param outerThresh A numeric value to threshold the distance from the local origin to the cell to classify cells as Outer cells.
#'@param wedges A numeric value for the number of wedge-shpaed bins for grouping the theta values.
#'@return A cell level datatable with new columns of XLocal, YLocal, RadialPostion, Theta, Outer, Density, Sparse and Perimeter values.
#'
#'@export
positionParms <- function(DT,outerThresh=.2, wedges=18, sparseThresh=.8){
  # count the number of cells within a euclidean distance from each cell

  #Add local cartesian and polar coordinates to each cell
  DT <- DT[,XLocal := X-median(X), by="Barcode,Well,Spot"]
  DT <- DT[,YLocal := Y-median(Y), by="Barcode,Well,Spot"]
  DT <- DT[,RadialPosition := sqrt(XLocal^2+YLocal^2)]
  DT <- DT[,Theta := calcTheta(XLocal,YLocal)]

  #Add the cell density
  #Average nuclear radius is 40 so touching nuclei are 80 apart
  #Set neighborhood as 4 nuclei radii
  DT <- DT[,Density:=spotCellDensities(.SD, radius=160)*10000,by="Barcode,Well,Spot"]
  DT <- DT[,Sparse := Density < sparseThresh]

  #Add a local wedge ID to each cell based on conversations with Michel Nederlof
  wedgeAngs <- 360/wedges
  DT <- DT[,Wedge:=ceiling(Theta/wedgeAngs)]

  #Define the perimeter cell if it exists in each wedge
  #Classify cells as outer if they have a radial position greater than a thresh
  DT <- DT[,OuterCell := labelOuterCells(RadialPosition, thresh=outerThresh),by="Barcode,Well,Spot"]

  #Require the cell not be in a sparse region
  denseOuterDT <- DT[!DT$Sparse  & DT$OuterCell]
  denseOuterDT <- denseOuterDT[,Perimeter := findPerimeterCell(.SD) ,by="Barcode,Well,Spot,Wedge"]
  setkey(DT,Barcode,Well,Object.ID)
  setkey(denseOuterDT,Barcode,Well,Object.ID)
  DT <- denseOuterDT[,list(Barcode,Well,Object.ID,Perimeter)][DT]
  DT$Perimeter[is.na(DT$Perimeter)] <- FALSE
  return(DT[,list(Barcode,Well,Spot,Object.ID,XLocal,YLocal,RadialPosition,Theta,Wedge,Density,Sparse,OuterCell,Perimeter)])
}


#'\code{calcTheta} Return the polar coordinate theta value
#'
#'@param x The x numeric value in cartesian coordinates.
#'@param y The y numeric value in cartesian coordinates.
#'@return The theta value of the X and y converted to polar coordinates.
#'
#'
#' Function from Roland on Stack Overflow
#'http://stackoverflow.com/questions/23018056/convert-cartesian-angles-to-polar-compass-cardinal-angles-in-r
calcTheta <- function(x,y) {
  z <- x + 1i * y
  res <- 90 - Arg(z) / pi * 180
  res %% 360
}

#'\code{findPerimeterCell} Determine the perimeter cell in wedge
#'
#' @param x A datatable or dataframe with a RadialPosition column
#' @return A logical vector the length of x with a TRUE value for the Perimeter cell
#'
#'
findPerimeterCell <- function(x){
  if(!nrow(x)==0){
    perimeterLogicals <- vector(length=nrow(x))
    perimeterLogicals[which.max(x$RadialPosition)] <- TRUE
  }
  return(perimeterLogicals)
}

#'\code{labelOuterCells} Determine if a cell is in Outer cell in the spot
#'
#' @param x A numeric vector of RadialPositions
#' @param thresh A numeric value to threshold x. The returned logical will be TRUE for cells with x values greater than thresh.
#' @return A logical vector the length of x with a TRUE value for the Outer cells
#'
#'
labelOuterCells <- function(x, thresh=.75){
  outerLogicals <- NULL
  if(!length(x)==0){
    outerLogicals <- x>quantile(x,probs = thresh, na.rm=TRUE)
  }
  return(outerLogicals)
}

