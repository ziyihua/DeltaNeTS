#' Solve the linear regression problem as formulated in DeltaNeTS with Lasso. Each gene is solved independently and sequentially.
#' @param LFC - log2FC expression data as matrix or dataframe, each row represents a gene and each column represents a sample
#' @param slope - slope matrix calculated by log2FC values at different time points within the same group of experiments
#' @param options - a list with fields 'alpha' and 'lambda'. Both alpha and lambda are sequence of values with alpha between 0 and 1 and lambda>=0.
#'                  See glmnet. The alpha field will be ignored because in this setting of Lasso, alpha is set to 1.
#' @param kfold - number of folds used in cross validation.
#' @param startG - index of gene to start with.
#' @param endG - index of the last gene (inclusive).
#' @return dlasso - a list with 3 fields
#'         @field A - A matrix
#'         @field P - Perturbation matrix
#'         @field Gindex - corresponding gene index

function(LFC,slope,options,kfold,startG,endG){
  n <- dim(LFC)[1]
  mtot <- dim(LFC)[2]

  X <- t(LFC)
  U <- diag(mtot)
  R <- t(slope)

  m1 <- dim(X)[1]
  m2 <- dim(R)[1]
  Z <- matrix(0,nrow = m2, ncol = dim(U)[2])

  n_df <- n
  grange <- seq(startG,endG)

  print("About to start DeltaNeTS-LASSO")

  norm_vec <- function(x) sqrt(sum(x^2))

  crossvalidate_glmnet_betaA_TS <- dget("crossvalidate_glmnet_betaA_TS.R")

  for (gi in seq(1,length(grange))) {
    gs <- grange[gi]
    cat(sprintf("DeltaNeTS-LASSO is running...(%4d/%d)\n",gi,length(grange)))

    if(sum(R[,gs])){
      r <- sqrt((norm_vec(X[,gs])^2/m1)/(norm_vec(R[,gs])^2/nnzero(R[,gs])))
    }else{
      r <- 1
    }

    Xts <- rbind(cbind(X,U),cbind(r*R,Z))
    yts <- c(X[,gs],r*R[,gs])

    options$exclude <- gs

    result <- crossvalidate_glmnet_betaA_TS(Xts,yts,n_df,options,kfold)

    if(gi==1){
      beta_res <- t(result$b_tmin)
    }else{
      beta_res <- cbind(beta_res,t(result$b_tmin))
    }

  }

  Amatrix <- t(beta_res[1:n,])
  Pmatrix <- t(beta_res[(n+1):dim(beta_res)[1],])

  print("Calculation is done.")
  dlasso <- list(A=Amatrix,P=Pmatrix,Gindex=grange)

  return(dlasso)

}
