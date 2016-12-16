#' 1-D linear interpolation
#' @param x - a monotonically increasing vector
#' @param y - a vector with length(x) entries or a matrix with length(x) rows
#' @param xi - a vector
#' @return yi - value of y at the points of xi

function(x,y,xi){
  
  if(dim(y)[1]!=length(x)){
    stop("Dimensions of x and y do not match.")
  }
  k <- order(xi)
  xxi <- xi[k]
  j <- order(c(x,xxi))
  dum <- c(x,xxi)[j]
  r <- c(1:length(j))[order(j)]
  r <- r[(length(x)+1):length(r)]-c(1:length(xxi))
  r[k] <- r
  r[which(xi==x[length(x)])] <- length(x)-1
  ind <- which((r>0)&(r<length(x)))
  yi <- matrix(NA,nrow = length(xxi),ncol = dim(y)[2])
  rind <- r[ind]
  u <- (xi[ind]-x[rind])/(x[rind+1]-x[rind])
  yi[ind,] <- y[rind,]+(y[rind+1,]-y[rind,])*t(tcrossprod(rep(1,dim(y)[2]),u))
  return(yi)
}