#' Load gene list file, gene regulatory network file (TF-gene interaction) and protein-protein interaction file
#' Please change the column based on the format of your network file accordingly
#' pdi_g: target genes in TF-gene interaction (gene regulatory network) file
#' pdi_tf: transcription factor in TF-gene interaction (gene regulatory network) file
#' ppi_tf: transcription factor in protein-protein interaction file
#' ppi_p: upstream regulator protein in protein-protein interaction file
#' You only need to change this part of the script, all other codes should be left unchanged
glist <- as.character(read.table("../example_data/GList_drug.txt")$V1)
pdi <- read.table('../example_data/13_lymphoma.txt')
ppi <- read.table('../example_data/tf_ppi_enrichr.txt')
pdi_g <- as.character(pdi$V2)
pdi_tf <- as.character(pdi$V1)
ppi_tf <- as.character(ppi$V1)
ppi_p <- as.character(ppi$V2)
#' Load Perturbation matrix
P <- read.table('../example_data/Pmat_drug.txt',sep = ',')


#' Remove target genes(pdi_g) that are not in the gene list
#' Add transcription factors to glist if some of their target genes are already in glist
#' Get index of pdi_g and pdi_tf in glist
pdi_g_index <- match(pdi_g,glist)
pdi_g <- pdi_g[!is.na(pdi_g_index)]
pdi_tf <- pdi_tf[!is.na(pdi_g_index)]
new_tf <- setdiff(unique(pdi_tf),glist)
glist <- c(glist,new_tf)
pdi_g_index <- match(pdi_g,glist)
pdi_tf_index <- match(pdi_tf,glist)

#' Remove transcription factors(ppi_tf) that are not in the gene list
#' Add upstream regulator proteins to glist if some of their interacting transcription factors are already in glist
#' Get index of ppi_tf and ppi_p in glist
ppi_tf_index <- match(ppi_tf,glist)
ppi_tf <- ppi_tf[!is.na(ppi_tf_index)]
ppi_p <- ppi_p[!is.na(ppi_tf_index)]
new_p <- setdiff(unique(ppi_p),glist)
glist <- c(glist,new_p)
ppi_tf_index <- match(ppi_tf,glist)
ppi_p_index <- match(ppi_p,glist)

#' Construct TF-gene and protein-TF interaction adjacency matrix 
#' To speed up matrix multiplication, the matrices use sparse notation
#' The final adjacency matrix by combining the aforementioned two matrices is calculated as (adj_ppi+I)*adj_pdi
library(slam)
adj_ppi <- matrix(0,nrow = length(glist),ncol = length(glist))
adj_ppi[cbind(ppi_p_index,ppi_tf_index)]=1
adj_ppi <- as.simple_triplet_matrix(t(adj_ppi+diag(length(glist))))
adj_pdi <- matrix(0,nrow = length(glist),ncol = length(glist))
adj_pdi[cbind(pdi_tf_index,pdi_g_index)]=1
adj_pdi <- as.simple_triplet_matrix(adj_pdi)
adj_tot <- crossprod_simple_triplet_matrix(adj_ppi,adj_pdi)

#' Calculate upstream perturbation matrix based on adjacency matrix
#' By averaging all contributions from its downstream effectors
#' Pmean is used as the final result
n <- length(glist)
m <- dim(P)[2]
Peval <- rbind(P,matrix(0,nrow = n-dim(P)[1], ncol = m))
Psum <- matrix(0,nrow = n, ncol = m)
Pmean <- matrix(0,nrow = n, ncol = m)
for (i in 1:n) {
  tgi <- which(adj_tot[i,]!=0)
  dn <- length(tgi)
  if(dn==0){
    dn <- 1
  }
  Psum[i,] <- apply(abs(Peval[tgi,]), 2, sum)
  Pmean[i,] <- Psum[i,]/dn
}