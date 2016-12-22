#' Convert a gene regulatory network with gene IDs to indices in the gene list and remove genes that are not in the gene list
#' @param tfnet - a tab delimited regulatory network file with gene names. The first column is the gene names of transcription
#'                factors and the second column is gene being regulated.
#' @param glist - gene ID file. The tfnet and glist file should use the same system for gene names. The genes in glist should
#'                have the same order as the expression data.
#' @param saveFile - TRUE or FALSE indicating whether to save the network file for the purpose of reuse. Default is FALSE. If 
#'                   not specified, the file will not be saved.
#' @return edge - a N*2 matrix where N is the number of TF-gene interactions. The first column is the index of TF and the second
#'                column is the index of the gene regulated by a TF. 

function(tfnet,glist,saveFile=FALSE){
  tfnet <- read.table(tfnet)
  glist <- as.character(read.table(glist)$V1)
  G1 <- as.character(tfnet$V1)
  G2 <- as.character(tfnet$V2)
  
  edge <- matrix(nrow = length(G1),ncol = 2)
  
  edge[,1] <- match(G1,glist)
  edge[,2] <- match(G2,glist)
  keep <- !is.na(edge[,1]) & !is.na(edge[,2])
  edge <- edge[keep,]
  keep <- edge[,1]!=edge[,2]
  edge <- edge[keep,]
  
  if(saveFile){
    write.table(edge,file = "../example_data/edge.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
  }
  
  return(edge)
}
