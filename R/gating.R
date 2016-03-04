#' Find a local minima appropriate for classiying a vector of numbers
#'
localMinima <- function(x, probs=c(.2,.8)){
  #browser()
  #Finds the local minima between the probs quantiles
  #x numeric vector
  #probs interval limits on where to search for the minima
  h <- hist(x,breaks=300, plot=FALSE)
  if(length(h$mids)<2) return(max(x))
  f <- approxfun(h$mids, h$counts)
  o <- optimise(f, interval=quantile(x, probs))
  if(length(o)>2) stop()
  return(o$minimum)
}

#' Gate a vector on a local minima
#' @export
gateOnlocalMinima <- function(x, ...){
  thresh <- MEMA:::localMinima(x, ...)
  cluster <- rep.int(1,times=length(x))
  cluster[x>thresh] <- 2
  return(cluster)
}
