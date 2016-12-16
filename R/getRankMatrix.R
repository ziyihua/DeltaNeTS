#' Calculate rank of genes based on perturbation matrix in each sample/experiment
#' @param Pmatrix - perturbation matrix, each row represents a gene and each column represents a sample
#' @param replist - a vector of length equal to #samples in Pmatrix, indicates which samples should be grouped together.
#'                  Default is NULL. If not specified, aggregation method will not be applied.
#' @param method - 'Mean' or 'Median'. Aggregation method to be applied to values of Pmatrix in each group of samples.
#' @return ave_result - a matrix of gene ranks, each row represents a gene and each column represents an aggregated sample

function(Pmatrix,replist=NULL,method='Median'){
  if(is.null(replist)|length(replist)==0){
    replist <- seq(1,dim(Pmatrix)[2])
  }else if(length(replist)!=dim(Pmatrix)[2]){
    stop("The length of replist should be the same as the number of samples (= the number of columns of LogFC data matrix.)")
  }
  
  rankMatrix <- matrix(1*dim(Pmatrix)[1],nrow = dim(Pmatrix)[1],ncol = max(replist)) # rank matrix

  for (i in seq(1,max(replist))) {
    ri <- which(replist==i)
    Psub <- Pmatrix[,ri]
    if(!is.null(dim(Psub))){
      if(method=='Median'){
        Pave <- apply(Pmatrix[,ri],1,median)
      }else if(method=='Mean'){
        Pave <- apply(Pmatrix[,ri],1,mean)
      }
    }else{
      Pave <- Psub
    }
    sorti <- order(abs(Pave),decreasing = TRUE)
    sortval <- abs(Pave)[sorti]
    rankMatrix[sorti[which(sortval>0)],i] <- seq(1:length(which(sortval>0)))
    if(i==1){
      Puni <- Pave
    }else{
      Puni <- cbind(Puni,Pave)
    }
  }
  
  ave_result <- list(P=Puni,R=rankMatrix)
  return(ave_result)
  
}