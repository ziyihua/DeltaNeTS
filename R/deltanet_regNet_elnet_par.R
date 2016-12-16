#' Solve the linear regression problem as formulated in DeltaNet with prior information on gene regulatory network and elastic net.
#' Each gene is solved independently and in parallel.
#' This function uses parallel computing.
#' @param LFC - log2FC expression data as matrix or dataframe, each row represents a gene and each column represents a sample
#' @param regNet - edge representation of a regulatory network. See 'CreateRegNet.R'.
#' @param options - a list with fields 'alpha' and 'lambda'. Both alpha and lambda are sequence of values with alpha between 0 and 1 and lambda>=0.
#'                  See glmnet.
#' @param kfold - number of folds used in cross validation.
#' @param startG - index of gene to start with.
#' @param endG - index of the last gene (inclusive).
#' @return dlasso - a list with 3 fields
#'         @field A - A matrix
#'         @field P - Perturbation matrix
#'         @field Gindex - corresponding gene index
#'         @field Gindex_solved - indices of genes that are solved (with RegNet it is possible that some genes are not solved)             
#'         @field alpha - alpha value that gives the smallest cv error for each gene
#'         @field lambda - lambda value that gives the smallest cv error for each gene. Note that for each gene, the smallest 
#'                         cv error is given by a pair of alpha and lambda value.

function(LFC,regNet,options,kfold,startG,endG){
  N <- dim(LFC)[1]
  M <- dim(LFC)[2]

  X <- cbind(t(LFC),diag(M))
  Y <- t(LFC)

  edge <- regNet
  tgi <- unique(edge[,2])
  tgi <- sort(tgi)

  grange <- intersect(tgi,startG:endG)

  print("about to start DeltaNet ElNet")

  #writeLines(c(""),"log.txt")
  #sink("log.txt",append = TRUE)

  comb <- function(x,...){
    mapply('cbind',x, ...,SIMPLIFY = FALSE)
  }

  res <- foreach (gi=1:length(grange),.combine='comb',.packages = c("glmnet","cvTools"))%dopar% {
    gs <- grange[gi]
    cat(sprintf("DeltaNet-ElNet is running...(%4d/%d)\n",gi,length(grange)))

    y <- Y[,gs]

    tfi <- edge[which(edge[,2]==gs),1]
    Gex <- setdiff(1:N,tfi)

    options_pr <- options
    options_pr$exclude <- union(gs,Gex)

    cv <- cvFolds(length(y),K=kfold)
    cvi <- cv$which[order(cv$subsets)]

    lambda <-  vector()
    cvmin <-  vector()
    beta <- matrix(data=0,nrow = dim(X)[2],ncol = length(options$alpha))

    for (alphaIdx in 1:length(options_pr$alpha)) {
      cvfit <- cv.glmnet(X,y,family='gaussian',alpha=options_pr$alpha[alphaIdx],lambda=options_pr$lambda,
                         exclude=options_pr$exclude,
                         intercept=0,foldid=cvi)
      lambda[alphaIdx] <- cvfit$lambda.min
      cvmin[alphaIdx] <- cvfit$cvm[cvfit$lambda==cvfit$lambda.min]
      coef.cvfit <- coef.cv.glmnet(cvfit,s='lambda.min')
      beta[coef.cvfit@i,alphaIdx] = coef.cvfit@x
    }

    minIdx <- max(which(cvmin==min(cvmin)))
    r <- beta[,minIdx]
    list(r,lambda[minIdx],options$alpha[minIdx])
  }
  #sink()

  beta_res <- res[[1]]
  Gindex <- which(startG:endG %in% grange)
  Amatrix <- matrix(data=0,nrow = endG-startG+1, ncol = N)
  Amat <- t(beta_res[1:N,])
  Amatrix[Gindex,] <- Amat
  Pmat <- t(beta_res[(N+1):dim(beta_res)[1],])
  Pmatrix <- LFC[startG:endG,]
  Pmatrix[Gindex,] <- Pmat

  print("Calculation is done.")
  dlasso <- list(A=Amatrix,P=Pmatrix,Gindex=(startG:endG),Gindex_solved=grange,alpha=res[[3]],lambda=res[[2]])

  return(dlasso)
}
