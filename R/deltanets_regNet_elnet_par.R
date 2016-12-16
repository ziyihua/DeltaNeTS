#' Solve the linear regression problem as formulated in DeltaNeTS with prior information on gene regulatory network and elastic net.
#' Each gene is solved independently and in parallel.
#' This function uses parallel computing.
#' @param LFC - log2FC expression data as matrix or dataframe, each row represents a gene and each column represents a sample
#' @param slope - slope matrix calculated by log2FC values at different time points within the same group of experiments
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

function(LFC,slope,regNet,options,kfold,startG,endG){
  n <- dim(LFC)[1]
  m_lfc <- dim(LFC)[2]
  m_slope <- dim(slope)[2]

  U <- diag(m_lfc)
  Z <- matrix(0,nrow = m_slope, ncol = m_lfc)

  n_df <- n
  delG <- which(apply(abs(slope), 1, sum)==0)
  index <- which(regNet[,1] %in% delG)
  edge <- regNet[-index,]
  tgi <- unique(edge[,2])
  tgi <- sort(tgi)

  grange <- intersect(tgi,startG:endG)

  norm_vec <- function(x) sqrt(sum(x^2))

  print("about to start DeltaNeTS ElNet")

  #writeLines(c(""),"log.txt")
  #sink("log.txt",append = TRUE)

  comb <- function(x,...){
    mapply('cbind',x, ...,SIMPLIFY = FALSE)
  }

  res <- foreach (gi=1:length(grange),.combine='comb',.packages = c("glmnet","cvTools"))%dopar% {
    gs <- grange[gi]
    cat(sprintf("DeltaNeTS-ElNet is running...(%4d/%d)\n",gi,length(grange)))

    if(sum(slope[gs,])){
      alpha <- sqrt((norm_vec(LFC[gs,])^2/m_lfc)/(norm_vec(slope[gs,])^2/nnzero(slope[gs,])))
    }else{
      alpha <- 1
    }

    Xin <- rbind(cbind(t(LFC),U),cbind(alpha*t(slope),Z))
    yin <- as.numeric(c(LFC[gs,],alpha*slope[gs,]))

    tfi <- edge[which(edge[,2]==gs),1]
    Gex <- setdiff(1:n,tfi)

    options_pr <- options
    options_pr$exclude <- union(gs,Gex)

    cv <- cvFolds(length(yin),K=kfold)
    cvi <- cv$which[order(cv$subsets)]

    lambda <-  vector()
    cvmin <-  vector()
    beta <- matrix(data=0,nrow = dim(Xin)[2],ncol = length(options$alpha))

    for (alphaIdx in 1:length(options_pr$alpha)) {
      cvfit <- cv.glmnet(Xin,yin,family='gaussian',alpha=options_pr$alpha[alphaIdx],lambda=options_pr$lambda,
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
  Amatrix <- matrix(data=0,nrow = endG-startG+1, ncol = n)
  Amat <- t(beta_res[1:n,])
  Amatrix[Gindex,] <- Amat
  Pmat <- t(beta_res[(n+1):dim(beta_res)[1],])
  Pmatrix <- LFC[startG:endG,]
  Pmatrix[Gindex,] <- Pmat
  print("Calculation is done.")
  dlasso <- list(A=Amatrix,P=Pmatrix,Gindex=(startG:endG),Gindex_solved=grange,alpha=res[[3]],lambda=res[[2]])

  return(dlasso)


}
