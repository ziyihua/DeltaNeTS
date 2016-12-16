#' Calculate slope at a time point with second order accurate finite difference approximations using 3 time points
#' @param tpoints - time points of each sample
#' @param f - samples at each time point where each column correspond to a time point
#' @param tsi - 1,2 or 3 indicating for which sample/time point the slope is calculated
#' @return dfdt - a vector of slope

findiff <- function(tpoints,f,tsi){
  tdiff <- tpoints - tpoints[tsi]
  if(length(which(tdiff>0))==0){
    deltsind <- max(which(tdiff<0))
    signdelts <- -1
  }else{
    deltsind <- min(which(tdiff>0))
    signdelts <- 1
  }
  delt <- tdiff/abs(tdiff[deltsind])

  Parind <- 1:length(tpoints)
  Parind <- Parind[-c(tsi,deltsind)]

  xmat <- delt[Parind]^2
  sol <- -1/xmat

  dfdt <- (-(1+sol)*f[,tsi]+f[,deltsind]+f[,Parind]*sol)/((signdelts+delt[Parind]*sol)*abs(tdiff[deltsind]))
  return(dfdt)
}
