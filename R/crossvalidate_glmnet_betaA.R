#' Solve for one gene as formulated in DeltaNet (Y=XB) with Lasso and cross validation. 
#' @param X - design matrix
#' @param y - vector of dependent variables
#' @param options - a list with fields 'alpha' and 'lambda'. Both alpha and lambda are sequence of values with alpha between 0 and 1 and lambda>=0.
#'                  See glmnet. The alpha field will be ignored because in this setting of Lasso, alpha is set to 1.
#' @param kfold - number of folds used in cross validation.
#' @return result - a list with 8 fields
#'         @field b_tmin - a vector of coefficients yielding the smallest cv error
#'         @field b_t1se - a vector of coefficients selected by 1-SE rule
#'         @field sopt_min - lambda value yielding the smallest cv error
#'         @field sopt_1se - lambda value selected by 1-SE rule
#'         @field cvm - mean of cv errors
#'         @field cvsd - standard deviation of cv errors
#'         @field b_series - a vector of coefficients given by glmnet
#'         @field tA - L1 norm of coefficients

function(X,y,options,kfold){
  cv <- cvFolds(length(y),K=kfold)
  cvi <- cv$which[order(cv$subsets)]
  srange <- seq(0,1,1e-3)


  cvk <- matrix(nrow = kfold, ncol = length(srange))
  for(k in seq(1,kfold)){

    foldi <- which(cvi==k)
    foldni <- which(cvi!=k)
    kfit <- glmnet(X[foldni,],y[foldni],family="gaussian",alpha=1,lambda = options$lambda,
                   exclude = options$exclude,intercept = 0)
    beta <- as(kfit$beta,"matrix")
    t <- apply(abs(beta[1:(dim(beta)[1]-length(y)),]),2,sum)
    s <- t/max(t)
    interp1q <- dget("interp1q.R")
    b_interpolated <- interp1q(s,t(beta),srange)

    Xtest <- X[foldi,]
    ytest <- y[foldi]

    residual <- matrix(rep(ytest,length(srange)),ncol = length(srange))-Xtest%*%t(b_interpolated)
    cvk[k,] <- apply(residual^2,2,sum)/length(ytest)
  }


  cvm <- apply(cvk,2,mean)
  cvsd <- apply(cvk,2,sd)
  if(!anyNA(cvm)){
    minval <- min(cvm)
    opti1 <- which.min(cvm)
    opti2 <- min(which(cvm<minval+cvsd[opti1]))
  }else{
    cvm_new <- cvm
    cvm_new[which(is.na(cvm_new))] <- Inf
    minval <- min(cvm_new)
    opti1 <- which.min(cvm_new)
    opti2 <- min(which(cvm_new<minval+cvsd[opti1]))
  }
  sopt1 <- srange[opti1]
  sopt2 <- srange[opti2]



  fit <- glmnet(X,y,family="gaussian",alpha=1,lambda = options$lambda,
                exclude = options$exclude,intercept = options$intr)
  beta <- as(fit$beta,"matrix")
  t <- apply(abs(beta[1:(dim(beta)[1]-length(y)),]),2,sum)
  s <- t/max(t)
  b_opt <- interp1q(s,t(beta),c(sopt1,sopt2))

  if(sopt1==1){
    b_opt[1,] <- beta[,dim(beta)[2]]
  }


  if(sopt2==1){
    b_opt[2,] <- beta[,dim(beta)[2]]
  }

  result.b_tmin <- t(b_opt[1,])
  result.b_t1se <- t(b_opt[2,])
  result.sopt_min <- sopt1
  result.sopt_1se <- sopt2
  result.cvm <- cvm
  result.cvsd <- cvsd
  result.b_series <- beta
  result.tA <- t
  result <- list(b_tmin=result.b_tmin,b_t1se=result.b_t1se,sopt_min=result.sopt_min,
                 sopt_1se=result.sopt_1se,cvm=result.cvm,cvsd=result.cvsd,
                 b_series=result.b_series,tA=result.tA)
  return(result)
}
