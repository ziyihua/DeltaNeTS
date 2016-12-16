#' Create edge representation of regulatory network
createRegNet <- dget("createRegNet.R")
regNet <- createRegNet("../example_data/32_lung_epithelium_lung_cancer.txt","../example_data/GList_full.txt")

#' Load files: log2fc, time points, group index, gene list and p-Value
lfc <- read.table("../example_data/LFC.txt")
tp <- read.table("../example_data/tpR.txt")$V1
group <- read.table("../example_data/grpR.txt")$V1
glist <- as.character(read.table("../example_data/GList_full.txt")$V1)
pval <- read.table("../example_data/pval.txt")

#' Sequence of lambda and alpha values used in grid search 
#' lambda and alpha are parameters to glmnet
opt.alpha <- c(10^seq(0,-2.5,length.out=10),0)
options <- list(alpha=opt.alpha)

#' Run DeltaNeTS and save the result
run_DeltaNeTS <- dget("run_DeltaNeTS.R")
result <- run_DeltaNeTS(lfc=lfc,tp=tp,group=group,regNet=regNet,pVal=pval,options=options,endG=5000,numCores = 10,parallel = TRUE)
save(result,file = "Calu3_1_5000")

#' Get top 100 genes for each experiments 
getRankMatrix <- dget("getRankMatrix.R")
rank <- getRankMatrix(P,sort(group),method = 'Mean')$R
topGenes <- matrix(nrow = 100,ncol = dim(rank)[2])
for (i in 1:dim(topGenes)[2]) {
  rankPerExp <- rank[,i]
  topIndx <- which(rankPerExp<=100)
  sorti <- order(rankPerExp[topIndx],decreasing = TRUE)
  topGenes[,i] <- glist[topIndx[sorti]]
}
write.table(topGenes[,15],file = "top.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
