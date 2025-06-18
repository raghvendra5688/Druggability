library(data.table)
library(ggplot2)

setwd("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/Druggability_Score/scripts/")

ppi_df <- fread("../../General_Data/9606.protein.links.full.v12.0.txt.gz",header=T)
ppi_df <- as.data.frame(ppi_df)

#Subset protein protein interaction with confidence score > 400
subset_ppi_df <- ppi_df[ppi_df$combined_score>=400,]

#Load the protein information
protein_info_df <- fread("../../General_Data/9606.protein.info.v12.0.txt.gz",header=T)
protein_info_df <- as.data.frame(protein_info_df)
protein_ensmbl_names <- protein_info_df$`#string_protein_id`

rev_ppi_df <- subset_ppi_df[subset_ppi_df$protein1 %in% protein_ensmbl_names,]
final_ppi_df <- rev_ppi_df[rev_ppi_df$protein2 %in% protein_ensmbl_names,]
final_ppi_df$gene1 <- NA
final_ppi_df$gene2 <- NA

for (i in 1:nrow(final_ppi_df))
{
  protein1 <- final_ppi_df$protein1[i]
  protein2 <- final_ppi_df$protein2[i]
  protein_name1 <- protein_info_df[protein_info_df$`#string_protein_id`==protein1,]$preferred_name
  protein_name2 <- protein_info_df[protein_info_df$`#string_protein_id`==protein2,]$preferred_name
  final_ppi_df$gene1[i] <- protein_name1
  final_ppi_df$gene2[i] <- protein_name2
}

#String Database: gene1, gene2, score
stringdb_edge_df <- final_ppi_df[,c("gene1","gene2","combined_score")]
write.table(stringdb_edge_df, file = "../Data/StringDB_interactions.csv",row.names=F, col.names=T, quote=F, sep=",")

#Subset Protein information
################################################################################
unique_genes <- unique(c(stringdb_edge_df$gene1,stringdb_edge_df$gene2))
final_protein_df <- protein_info_df[protein_info_df$preferred_name%in%unique_genes,]
write.table(final_protein_df[,c(1:3)], file = "../Data/StringDB_protein_info.csv",row.names=F, col.names=T, quote=F, sep=" | ")
