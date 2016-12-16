#' Choose the correct setting of DeltaNeTS to run and construct a time slope matrix based on your input arguments
#' @param lfc - log2FC expression data as matrix or dataframe, each row represents a gene and each column represents a sample
#' @param tp - timepoint of each sample. Default is NULL. If not specified, method 'DeltaNet' will be chosen.
#' @param group - group index of each sample. The same sample/experiment at different time points should have the same group 
#'                index. Repitition experiment should have different group index. Together with time point, this information 
#'                is needed to calculate time slope matrix. Default is NULL. If not specified, method 'DeltaNet' will be chosen.
#' @param regNet - edge representation of a regulatory network. See 'CreateRegNet.R'. Default is NULL. If not specified, 
#'                 DeltaNet/DeltaNeTS with Lasso will be chosen.
#' @param pVal - a p-Value matrix or dateframe representing the statistical significance of each log2FC value. It should have the same dimension
#'               as lfc. It is used to set some entries of the time slope matrix to zero. Default is null. If not specified, all slopes
#'               are used in the time series setting.
#' @param pValSig - a threshold between 0 and 1. If pVal is available, for each gene and all samples in a group, the time slope values are 
#'                  set to zero, if any of them has p-Value greater than pValSig. Default is 0.5.
#' @param options - a list with fields 'alpha' and 'lambda'. Both alpha and lambda are sequence of values with alpha between 0 and 1 and lambda>=0.
#'                  See glmnet. Default is NULL. If not specified, or only one field is specified, a default sequence is provided.
#'                  Note that if you don't provide regNet, DeltaNeTS/DeltaNet will solve Lasso and no matter what value you specify for alpha, it is
#'                  always set to 1.
#' @param parallel - TRUE or FALSE indicating whether to do parallel computing. Default is FALSE.
#' @param numCores - number of cores to use for parallel computing. Default is 4. This only matters if you set parallel to TRUE.
#' @param kfold - number of folds used in cross validation. Default is 10.
#' @param startG - index of gene to start with. Default is NULL. If not specified, it is set to 1.
#' @param endG - index of the last gene (inclusive). Default is NULL. If not specified, all genes from the startG to the last one will be solve.
#'               Since DeltaNeTS/DeltaNet solve for each gene independently, with startG and endG you can specify a range of gene to solve.                  
#' @return dlasso - a list with 3 or 6 fields depending on whether regNet is provided.
#'         @field A - A matrix
#'         @field P - Perturbation matrix. Note that if the method is DeltaNet, order of the samples in P is same as in lfc.
#'                    If the method is DeltaNeTS, the samples are sorted first by an increasing order of group index and then time points.
#'         @field Gindex - corresponding gene index
#'         with regNet there are three other fields:
#'         @field Gindex_solved - indices of genes that are solved (with RegNet it is possible that some genes are not solved)             
#'         @field alpha - alpha value that gives the smallest cv error for each gene
#'         @field lambda - lambda value that gives the smallest cv error for each gene. Note that for each gene, the smallest 
#'                         cv error is given by a pair of alpha and lambda value.

function(lfc,tp=NULL,group=NULL,regNet=NULL,pVal=NULL,pValSig=0.5,options=NULL,parallel=FALSE,numCores=4,kfold=10,startG=NULL,endG=NULL){

  #configure default parameters if not given by user
  if(is.null(options)){
    options <- list(alpha=c(10^seq(0,-2,length.out=20),0),lambda=10^seq(-0.1,-5,length.out = 100))
  }else{
    if(is.null(options$lambda)){
      options$lambda=10^seq(-0.1,-5,length.out = 100)
    }else if(length(options$lambda)<2){
      stop("More than one value of lambda is required for glmnet with CV.")
    }
    if(is.null(options$alpha)){
      options$alpha=c(10^seq(0,-2,length.out=20),0)
    }else if(any(options$alpha<0||any(options$alpha>1))){
      stop("Alpha should be in the range 0 and 1.")
    }
  }

  if(kfold%%1!=0|kfold<5){
    stop("kfold should be a natural number more than 5. Default value is 10.")
  }

  if(pValSig<=0|pValSig>1){
    stop("Threshold for pValue should be between 0 and 1.")
  }

  if(is.null(startG)){
    startG=1
  }else if(startG<1){
    stop("Gene index should be greater or equal than 1.")
  }
  if(is.null(endG)){
    endG=dim(lfc)[1]
  }else if(endG>dim(lfc)[1]){
    stop("Gene index exceeds dimension of the expression data.")
  }
  if(startG>endG){
    stop("startG should be less or equal to endG.")
  }


  if(is.null(group)||is.null(tp)){
    method <- "DeltaNet"
  }else{
    #reorder group
    if(length(group)!=dim(lfc)[2]||length(tp)!=dim(lfc)[2]){
      stop("Dimension of group or timepoint does not match number of samples of the expression data.")
    }
    grp_ord <- order(group)
    group <- group[grp_ord]
    tp <- tp[grp_ord]
    lfc <- lfc[,grp_ord]
    if(!is.null(pVal)){
      pVal <- pVal[,grp_ord]
    }

    #reorder time points within each group
    method <- "DeltaNet"
    grp_unique <- unique(group)

    if(length(grp_unique)!=length(group)){
      slopeMat <- vector(mode="numeric",length = 0L)
      findiff <- dget("findiff.R")
      for (i in grp_unique) {
        grp_index <- which(group==i)
        if(length(unique(tp[grp_index]))!=1){
          if(length(unique(tp[grp_index]))!=length(tp[grp_index])){
            stop("Replicated samples of some time points in one group. Samples in each group shoud have distinct time points or the same time point.")
          }
          tp_order <- order(tp[grp_index])
          grp_index_order <- grp_index[tp_order]
          tp[grp_index] <- tp[grp_index_order]
          lfc[,grp_index] <- lfc[,grp_index_order]
          if(!is.null(pVal)){
            pVal[,grp_index] <- pVal[,grp_index_order]
          }
          method <- "DeltaNeTS"
          current_tp <- tp[grp_index]
          current_lfc <- lfc[,grp_index]
          if(!is.null(pVal)){
            current_pVal <- apply(pVal[,grp_index],1,max)
          }

          l <- length(current_tp)
          if(l>2){
            for (j in 1:l) {
              if(j==1){
                slop_sp <- findiff(current_tp[1:3],current_lfc[,1:3],1)
              }else if(j==l){
                slop_sp <- findiff(current_tp[(l-2):l],current_lfc[,(l-2):l],3)
              }else{
                slop_sp <- findiff(current_tp[(j-1):(j+1)],current_lfc[,(j-1):(j+1)],2)
              }
              if(!is.null(pVal)){
                slop_sp[which(current_pVal>pValSig)]=0
              }
              slopeMat <- cbind(slopeMat,slop_sp)
            }
          }else{
            slop_sp <- (current_lfc[,2]-current_lfc[,1])/(current_tp[2]-current_tp[1])
            if(!is.null(pVal)){
              slop_sp[which(current_pVal>pValSig)]=0
            }
            slopeMat <- cbind(slopeMat,slop_sp)
          }
        }
      }
    }
  }
  library(glmnet)
  library(cvTools)
  library(methods)
  library(Matrix)

  if(parallel){
    library(doParallel)
    library(foreach)
    cl <- makeCluster(numCores,outfile='')
    registerDoParallel(cl)
  }
  ptm <- proc.time()
  if(method=='DeltaNet'){
    if(parallel){
      if(is.null(regNet)){
        deltanet_lasso_par <- dget("deltanet_lasso_par.R")
        dlasso <- deltanet_lasso_par(lfc,options,kfold,startG,endG)
      }else{
        deltanet_regNet_elnet_par <- dget("deltanet_regNet_elnet_par.R")
        dlasso <- deltanet_regNet_elnet_par(lfc,regNet,options,kfold,startG,endG)
      }
    }else{
      if(is.null(regNet)){
        deltanet_lasso <- dget("deltanet_lasso.R")
        dlasso <- deltanet_lasso(lfc,options,kfold,startG,endG)
      }else{
        deltanet_regNet_elnet <- dget("deltanet_regNet_elnet.R")
        dlasso <- deltanet_regNet_elnet(lfc,regNet,options,kfold,startG,endG)
      }
    }
  }else{
    if(parallel){
      if(is.null(regNet)){
        deltanets_par <- dget("deltanets_par.R")
        dlasso <- deltanets_par(lfc,slopeMat,options,kfold,startG,endG)
      }else{
        deltanets_regNet_par <- dget("deltanets_regNet_elnet_par.R")
        dlasso <- deltanets_regNet_par(lfc,slopeMat,regNet,options,kfold,startG,endG)
      }
    }else{
      if(is.null(regNet)){
        deltanets <- dget("deltanets.R")
        dlasso <- deltanets(lfc,slopeMat,options,kfold,startG,endG)
      }else{
        deltanets_regNet <- dget("deltanets_regNet_elnet.R")
        dlasso <- deltanets_regNet(lfc,slopeMat,regNet,options,kfold,startG,endG)
      }
    }
  }

  if(parallel){
    stopCluster(cl)
  }
  print(proc.time()-ptm)
  return(dlasso)

}
