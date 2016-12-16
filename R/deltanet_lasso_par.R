#' Solve the linear regression problem as formulated in DeltaNet with Lasso. Each gene is solved independently and in parallel.
#' This function uses parallel computing.
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

  grange <- seq(startG,endG)
  print("about to start DeltaNet Lasso")

  crossvalidate_glmnet_betaA <- dget("crossvalidate_glmnet_betaA.R")
  interp1q <- dget("interp1q.R")
  #writeLines(c(""),"log.txt")
  #sink("log.txt",append = TRUE)
  beta_res <- foreach (gi=1:length(grange),.combine=cbind,.packages = c("glmnet","cvTools"),.export="interp1q")%dopar% {
    gs <- grange[gi]
    cat(sprintf("DeltaNet-LASSO is running...(%4d/%d)\n",gi,length(grange)))

    y <- Y[,gs]
    options_pr <- options
    options_pr$exclude <- gs

    result <- crossvalidate_glmnet_betaA(X,y,options_pr,kfold)

    r <- t(result$b_tmin)
  }
  #sink()

  Amatrix <- t(beta_res[1:N,])
  Pmatrix <- t(beta_res[(N+1):dim(beta_res)[1],])

  print("Calculation is done.")
  dlasso <- list(A=Amatrix,P=Pmatrix,Gindex=grange)

  return(dlasso)
}
