#Functions to support MEMAs printed in 8 well plates

#' Calculate the neighborhood density around each cell
#'@param spot A datatable with X and Y columns
#'@param radius The radial distance that defines the neighborhood around the cell
#'@return A numeric vector of length nrow(spot) with the cell density values
#'
#'@export
spotCellDensitiesScanR <- function(spot,radius=(max(spot$X)-min(spot$X))/5) {
  distMatrix <- as.matrix(dist(spot[,list(X,Y)]))
  count <- apply(distMatrix, 2, function(x){sum(x <= radius) - 1})
  cellDensity <- count/(pi*radius^2)
  return(cellDensity)
}

#' Calculate the neighborhood density around each cell
#'@param spot A datatable with X and Y columns
#'@param radius The radial distance that defines the neighborhood around the cell
#'@return A numeric vector of length nrow(spot) with the cell density values
#'
#'@export
spotCellDensities<- function (spot, radius = (max(spot$Cells_Location_Center_X) - min(spot$Cells_Location_Center_X))/5) 
{
  distMatrix <- as.matrix(dist(spot[, list(Cells_Location_Center_X, Cells_Location_Center_Y)]))
  count <- apply(distMatrix, 2, function(x) {
    sum(x <= radius) - 1
  })
  cellDensity <- count/(pi * radius^2)
  return(cellDensity)
}

#'Deriving Position-based Parameters
#'
#' Calculate the position-based derived parameters for cell-level data to aid population heterogenaity studies.
#'@param DT A cell level data.table with plate-level x,y coordinates
#'@param densityRadius The radius of the circle around each nuclei defining its neighborhood.
#'@param outerThresh A quantile value between 0 and 1 used to threshold RadialPosition. The returned logical in the Outer column will be TRUE for cells with RadialPosition values in quantiles greater than outerThresh.
#'@param wedges A numeric value for the number of wedge-shaped bins for grouping the theta values.
#'@param sparseThresh A numeric value used to threshold the Density values. The returned logical in the Sparse column will be TRUE for cells with Density values greater than sparseThresh.
#'@return A cell level data.table of length nrow(DT) with columns of XLocal (numeric), YLocal (numeric), RadialPostion (numeric), Theta (numeric), Outer (logical), Density (numeric), Sparse (logical) and Perimeter (logical) values. XLocal and YLocal are cartesian coordinates of each nuclei with the origin at the median X and median Y of each spot. RadialPosition and Theta are polar coordinates for XLocal and YLocal. Outer is a logical for cells in quantiles greater that outerThresh. Density is a scaled value that represents the number of nuclei centers within a densityRadius of each cell. Sparse is a logical for whether each cell is in a neighborhood with a Density value less than sparseThresh. Perimeter is a logical for the cell meeting the following criteria: Outer, furthest from the origin in their wedge and not Sparse.
#'
#' Theta values function is from Roland on Stack Overflow
#'http://stackoverflow.com/questions/23018056/convert-cartesian-angles-to-polar-compass-cardinal-angles-in-r
#'
#'@export
positionParms <- function(DT,densityRadius = 160, outerThresh=.2, wedges=18, sparseThresh=.8){
  # count the number of cells within a euclidean distance from each cell

  lDT <- copy(DT)
  #Add local cartesian and polar coordinates to each cell
  lDT <- lDT[,XLocal := X-median(X), by="Barcode,Well,Spot"]
  lDT <- lDT[,YLocal := Y-median(Y), by="Barcode,Well,Spot"]
  lDT <- lDT[,RadialPosition := sqrt(XLocal^2+YLocal^2)]
  lDT <- lDT[,Theta := calcTheta(XLocal,YLocal)]

  #Add the cell density
  #Average nuclear radius is 40 so touching nuclei are 80 apart
  #Set neighborhood as 4 nuclei radii
  lDT <- lDT[,Density:=spotCellDensities(.SD, radius=densityRadius)*10000,by="Barcode,Well,Spot"]
  lDT <- lDT[,Sparse := as.logical(Density < sparseThresh)]

  #Add a local wedge ID to each cell based on conversations with Michel Nederlof
  wedgeAngs <- 360/wedges
  lDT <- lDT[,Wedge:=ceiling(Theta/wedgeAngs)]

  #Define the perimeter cell if it exists in each wedge
  #Classify cells as outer if they have a radial position greater than a thresh
  lDT <- lDT[,OuterCell := labelOuterCells(RadialPosition, thresh=outerThresh),by="Barcode,Well,Spot"]

  #Require the cell not be in a sparse region
  denseOuterDT <- lDT[!lDT$Sparse  & lDT$OuterCell]
  denseOuterDT <- denseOuterDT[,Perimeter := findPerimeterCell(.SD) ,by="Barcode,Well,Spot,Wedge"]
  setkey(lDT,Barcode,Well,Spot,ObjectID)
  setkey(denseOuterDT,Barcode,Well,Spot,ObjectID)
  lDT <- denseOuterDT[,list(Barcode,Well,Spot,ObjectID,Perimeter)][lDT]
  lDT$Perimeter[is.na(lDT$Perimeter)] <- FALSE
  return(lDT[,list(Barcode,Well,Spot,ObjectID,XLocal,YLocal,RadialPosition,Theta,Wedge,Density,Sparse,OuterCell,Perimeter)])
}

#'Calculate the polar coordinate theta value
#'
#'Return the polar coordinate theta value that starts at 0 along the positive Y axis and increases to 360 in a clockwise direction.
#'
#'@param x The x numeric value in cartesian coordinates.
#'@param y The y numeric value in cartesian coordinates.
#'@return The theta value of the X and y converted to polar coordinates.
#
#'@export
# Function from Roland on Stack Overflow
#http://stackoverflow.com/questions/23018056/convert-cartesian-angles-to-polar-compass-cardinal-angles-in-r
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
#' @export
findPerimeterCell <- function(x){
  if(!nrow(x)==0){
    perimeterLogicals <- vector(length=nrow(x))
    perimeterLogicals[which.max(x$RadialPosition)] <- TRUE
  }
  return(perimeterLogicals)
}

#'Classify Outer cells
#'
#'\code{labelOuterCells} Determine if a cell is an Outer cell in the spot
#'
#' @param x A numeric vector of RadialPositions
#' @param thresh A quantile value between 0 and 1 used to threshold x. The returned logical will be TRUE for cells with x values in quantiles greater than thresh.
#' @return A logical vector the length of x with a TRUE value for the Outer cells
#'
#' @export
labelOuterCells <- function(x, thresh=.75){
  outerLogicals <- NULL
  if(!length(x)==0){
    outerLogicals <- x>quantile(x,probs = thresh, na.rm=TRUE)
  }
  return(outerLogicals)
}

#' Cluster using kmeans
#'
#' \code{kmeansClusterValue} is a wrapper function for perfoming kmeans clustering
#' @param x A numeric vector to be clustered
#' @param centers The number of centers or clusters to find.
#' @return The cluster assignments for x using the base kmeans command.
#'
#' @export
kmeansClusterValue <- function (x, centers = 2)
{
  #browser()
  x <- data.frame(x)
  xkmeans <- kmeans(x, centers = centers)
  if(centers==2){
    if(xkmeans$centers[1] > xkmeans$centers[2]){
      tmp <- xkmeans$cluster == 1
      xkmeans$cluster[xkmeans$cluster == 2] <- 1L
      xkmeans$cluster[tmp] <- 2L
    }
  }
  return(xkmeans[["cluster"]])
}

#' Cluster using kmeans
#'
#' \code{kmeansDNAClusterValue} is a wrapper function for perfoming kmeans clustering
#' @param x A numeric vector to be clustered
#' @param centers The number of centers or clusters to find.
#' @return The cluster assignments for x using the base kmeans command.
#'
#' @export
kmeansDNACluster <- function (x, centers = 2) 
{
  #browser()
  x <- data.frame(x)
  xkmeans <- kmeans(x, centers = centers)
  #Swap cluster IDs to make sure cluster 2 has higher values
  if(centers==2){
    if(xkmeans$centers[1] > xkmeans$centers[2]){
      tmp <- xkmeans$cluster == 1
      xkmeans$cluster[xkmeans$cluster == 2] <- 1L
      xkmeans$cluster[tmp] <- 2L
    }
  }
  return(xkmeans$cluster)
}

#'kMeans cluster an intensity value
#'
#'Returns cluster assignments for intensity values based on a control ligand
#'@param x A dataframe or data.table
#'@param value A character vector of one column name in x that conatins non-negative intensity values
#'@param ctrlLigand A character vector of one Ligand name in x
#'
#'@return Cluster assignments of length nrow(x)for the value column based on a kmeans threshold calculated from the ligand values. The values are 1 for intensities above the threshold and 0 for values equal to or less than the trheshold.
#'
#'@export
kmeansCluster <- function(x,value,ctrlLigand="HighSerum"){
  ctrlClusters <- kmeansClusterValue(log2(1+x[[value]][grepl(ctrlLigand,x$Ligand)]))
  ctrlPositiveThresh <- min(x[[value]][grepl(ctrlLigand,x$Ligand)][ctrlClusters==2])
  clusters <- rep.int(0,nrow(x))
  clusters[x[[value]]>ctrlPositiveThresh] <- 1
  return(clusters)
}

#' Count the number of neigbors around a nuclei
#' 
#'@param spot A data.table with \code{spot$Nuclei_CP_AreaShape_Center_X} and 
#'\code{spot$Nuclei_CP_AreaShape_Center_Y} values for each cell. 
#'@param radius The radius around the center of each cell used to count neighbor cells.
#'@export

cellNeighbors<- function (spot, radius = (max(spot$Nuclei_CP_AreaShape_Center_X) - min(spot$Nuclei_CP_AreaShape_Center_X))/5) 
{
  distMatrix <- as.matrix(dist(spot[, list(Nuclei_CP_AreaShape_Center_X, Nuclei_CP_AreaShape_Center_Y)]))
  count <- apply(distMatrix, 2, function(x) {
    sum(x <= radius) - 1
  })
  return(count)
}

#' Calculate the proportion of cells that are classfied as having 2N DNA
#' 
#' @param x numeric vector of cycle states with values of 1 for 2n and 2 for 4N
#' @return proportion of cells in x that are in 2N
#' @export
calc2NProportion <- function(x){
  if(!length(x)) stop("Calculating 2N/4N proportion on an empty group")
  if(sum(grepl("[^12]",x)))stop("Invalid cycle state passed to calc2NProportion")
  if(sum(x==1)) {
    proportion2N <- sum(x==1)/length(x)
  } else proportion2N <- 0
  return(proportion2N)
}

#'
#'@export
#Calculate the proportions in the gates
calcProportion <- function(x){
  sum(x)/length(x)
}