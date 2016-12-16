#' Solve the linear regression problem as formulated in DeltaNet with Lasso. Each gene is solved independently and sequentially.
#' @param LFC - log2FC expression data as matrix or dataframe, each row represents a gene and each column represents a sample
#' @param options - a list with fields 'alpha' and 'lambda'. Both alpha and lambda are sequence of values with alpha between 0 and 1 and lambda>=0.
#'                  See glmnet. The alpha field will be ignored because in this setting of Lasso, alpha is set to 1.
#' @param kfold - number of folds used in cross validation.
#' @param startG - index of gene to start with.
#' @param endG - index of the last gene (inclusive).
#' @return dlasso - a list with 3 fields
#'         @field A - A matrix
#'         @field P - Perturbation matrix
#'         @field Gindex - corresponding gene index

function(LFC,options,kfold,startG,endG){
  N <- dim(LFC)[1]
  M <- dim(LFC)[2]

  X <- cbind(t(LFC),diag(M))
  Y <- t(LFC)

  crossvalidate_glmnet_betaA <- dget("crossvalidate_glmnet_betaA.R")
  
  print("about to start DeltaNet Lasso")

  grange <- seq(startG,endG)
  cat(sprintf("DeltaNet-LASSO is running...(0/%d)\n",length(grange)))

  for (gi in seq(1,length(grange))) {
    gs <- grange[gi]
    cat(sprintf("DeltaNet-LASSO is running...(%4d/%d)\n",gi,length(grange)))

    y <- Y[,gs]
    options$exclude <- gs
    result <- crossvalidate_glmnet_betaA(X,y,options,kfold)
    if(gi==1){
      beta_res <- t(result$b_tmin)
    }else{
      beta_res <- cbind(beta_res,t(result$b_tmin))
    }
  }


  Amatrix <- t(beta_res[1:N,])
  Pmatrix <- t(beta_res[(N+1):dim(beta_res)[1],])

  print("Calculation is done.")
  dlasso <- list(A=Amatrix,P=Pmatrix,Gindex=grange)

  return(dlasso)

}
