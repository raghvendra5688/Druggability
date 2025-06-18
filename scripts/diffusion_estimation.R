library(data.table)
library(ggplot2)
library(igraph)
library(devtools)
library(dnet)
library(R.utils)
library(doParallel)
registerDoParallel(20)

setwd("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/Druggability_Score/")

#Load the drug protein interaction edgelist from the Data folder
################################################################################
drug_protein_targets_df <- fread("Data/drug_protein_interactions.csv")
drug_protein_targets_df <- as.data.frame(drug_protein_targets_df)
unique_drugs <- drug_protein_targets_df$drug_id
N_drugs <- length(unique_drugs)+1

#Load the string ppi matrix as an adjacency matrix
################################################################################
string_ppi_df <- read.table("Data/subset_StringDB_interactions.csv",header=T,sep=",")
g_ppi <- graph_from_edgelist(as.matrix(string_ppi_df[,c(1:2)]), directed = T)
E(g_ppi)$weight <- string_ppi_df$combined_score/1000
N_ppi <- length(V(g_ppi))
all_genes <- V(g_ppi)$name

#Get the data frame of each drug and the affinity to all proteins from their target genes 
################################################################################
p0_matrix <- matrix(0, nrow=N_ppi, ncol=N_drugs)
rownames(p0_matrix) <- V(g_ppi)$name
all_drug_targets <- NULL
for (i in 1:(N_drugs))
{
  if (i<N_drugs)
  {
    drug_id <- unique_drugs[i]
    targets <- drug_protein_targets_df[drug_protein_targets_df$drug_id==drug_id,]$collated_gene_name
    targets <- unlist(strsplit(targets,split=" ; "))
    all_drug_targets <- union(all_drug_targets, targets)
    p0_vec <- rep(0, N_ppi)
    names(p0_vec) <- all_genes
    p0_vec[names(p0_vec) %in% targets] <- 1
    p0_matrix[,i] <- p0_vec
  }else{
    p0_vec <- rep(0,N_ppi)
    names(p0_vec) <- all_genes
    p0_vec[names(p0_vec) %in% all_drug_targets] <- 1
    p0_matrix[,i] <- p0_vec
  }
}

#This is the random walk with restart distance of all genes from target genes of a drug
pinf_matrix <- sqrt(dRWR(g = g_ppi, normalise = "laplacian", setSeeds = p0_matrix, restart = 0.25, multicores = 10))
rownames(pinf_matrix) <- V(g_ppi)$name
colnames(pinf_matrix) <- c(unique_drugs,"all_drugs")
pinf_matrix <- as.matrix(pinf_matrix)   #Obtained the local +global propagation matrix i.e. pinf_matrix
write.table(pinf_matrix, file="Results/drug_target_propagation_scores.csv",row.names=T, col.names=T, quote=F, sep="\t")

system(paste0("zip Results/Propagation_Scores.zip Results/drug_target_propagation_scores.csv"))
system(paste0("rm Results/drug_target_propagation_scores.csv"))
system(paste0("split -b 5M Results/Propagation_Scores.zip; rm Results/Propagation_Scores.zip"))
system(paste0("mv x* Results/Propagation_Scores/"))





