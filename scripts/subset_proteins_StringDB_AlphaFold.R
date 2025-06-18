library(data.table)
library(ggplot2)

setwd("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/Druggability_Score/scripts/")

#Load the stringdb protein dataset
###############################################################################
stringdb_protein_info_df <- fread("../Data/StringDB_protein_info.csv")
stringdb_protein_info_df <- as.data.frame(stringdb_protein_info_df)
gene_names_stringdb <- stringdb_protein_info_df$preferred_name

#Load the proteins with structure in AlphaFold repository
###############################################################################
alphafold_df <- fread("../Data/Curated_AlphaFold_UniProtKB.csv",header=T)
alphafold_df <- as.data.frame(alphafold_df)

#Get gene names and match with Alphafold subset
###############################################################################
gene_names_alphafold <- unlist(lapply(strsplit(alphafold_df$Gene_Name,split="_"),`[[`,1))
alphafold_df$hgnc_symbol <- gene_names_alphafold
intersect_gene_names <- intersect(gene_names_alphafold, gene_names_stringdb)    #There is low intersection suggesting that StringDB gene names are not matching with Alphafold gene names 

#Need a strategy to map StringDB genes to AlphaFold gene names
################################################################################
#Get all the human uniport ids
###############################################################################
all_uniprot_human_df <- fread("../../General_Data/uniprotkb_Human_2024_06_03.tsv")
all_uniprot_human_df <- as.data.frame(all_uniprot_human_df)

subset_uniprot_human_df <- all_uniprot_human_df[grep("Homo sapiens",all_uniprot_human_df$Organism),]
subset_uniprot_human_with_gene_names_df <- subset_uniprot_human_df[subset_uniprot_human_df$`Gene Names`!="",]

#For each StringDB gene name find the corresponding set of Uniprot ids
################################################################################
mapping_df <- NULL
for (i in 1:length(gene_names_stringdb))
{
  gene_name <- gene_names_stringdb[i]
  uniprot_ids <- paste(subset_uniprot_human_with_gene_names_df[grep(gene_name,subset_uniprot_human_with_gene_names_df$`Gene Names`),]$Entry,collapse=" ; ")
  temp <- cbind(gene_name, uniprot_ids)
  mapping_df <- rbind(mapping_df, temp)
}
mapping_df <- as.data.frame(mapping_df)
colnames(mapping_df) <- c("gene_name_stringdb","uniprot_ids")

#Get the subset of genes with at least one uniprot id
################################################################################
subset_mapping_df <- mapping_df[which(sapply(strsplit(mapping_df$uniprot_ids,split=" ; "),length)>0),]

#Look at all Uniprot ids for human in Alphafold database and find matching gene names in the mapping dataframe
################################################################################
matched_gene_alphafold_uniprot_df <- NULL
for (i in 1:nrow(alphafold_df))
{
  alphafold_uniprot_id <- alphafold_df$UniProtId[i]
  protein_length <- subset_uniprot_human_with_gene_names_df[which(subset_uniprot_human_with_gene_names_df==alphafold_uniprot_id),]$Length
  temp_ids <- grep(alphafold_uniprot_id, subset_mapping_df$uniprot_ids)
  if (length(temp_ids)>0)
  {
    gene_names <- subset_mapping_df$gene_name_stringdb[temp_ids]
    temp <- cbind(gene_names,rep(alphafold_uniprot_id,length(gene_names)),rep(protein_length,length(gene_names)))
    matched_gene_alphafold_uniprot_df <- rbind(matched_gene_alphafold_uniprot_df, temp)
  }
}
matched_gene_alphafold_uniprot_df <- as.data.frame(matched_gene_alphafold_uniprot_df)
colnames(matched_gene_alphafold_uniprot_df) <- c("StringDB_Gene_Name","AlphaFold_Uniprot_Id","Uniprot_Protein_Length")
matched_gene_alphafold_uniprot_df$Uniprot_Protein_Length <- as.numeric(as.vector(matched_gene_alphafold_uniprot_df$Uniprot_Protein_Length))
matched_gene_alphafold_uniprot_df <- matched_gene_alphafold_uniprot_df[order(matched_gene_alphafold_uniprot_df$StringDB_Gene_Name,matched_gene_alphafold_uniprot_df$Uniprot_Protein_Length,decreasing=T),]
write.table(matched_gene_alphafold_uniprot_df, "../Data/StringDB_genes_AlphaFold_Mapping.csv",row.names=F, col.names=T, quote=F, sep="\t")

#Filter the Alphafold proteins to those which map to atleast one gene name
################################################################################
unique_alphafold_proteins <- unique(matched_gene_alphafold_uniprot_df$AlphaFold_Uniprot_Id)
filtered_alphafold_df <- alphafold_df[alphafold_df$UniProtId %in% unique_alphafold_proteins,c(1,2,3,5)]
write.table(filtered_alphafold_df, file="../Data/Filtered_AlphaFold_UniProtKB.csv",row.names=F, col.names=T, quote=F, sep=",")

#Find the subset of human genes with at least one uniprot id
subset_stringdb_protein_info_df <- stringdb_protein_info_df[stringdb_protein_info_df$preferred_name %in% matched_gene_alphafold_uniprot_df$StringDB_Gene_Name,]
write.table(subset_stringdb_protein_info_df,file = "../Data/subset_StringDB_protein_info.csv", row.names=F, col.names=T, quote=F, sep=" | ")

#Read the stringdb interactions and keep it to the subset of genes with one uniprot id
stringdb_ppi_df <- fread("../Data/StringDB_interactions.csv")
stringdb_ppi_df <- as.data.frame(stringdb_ppi_df)
subset_stringdb_ppi_df <- stringdb_ppi_df[stringdb_ppi_df$gene1 %in% matched_gene_alphafold_uniprot_df$StringDB_Gene_Name & stringdb_ppi_df$gene2 %in% matched_gene_alphafold_uniprot_df$StringDB_Gene_Name,]
write.table(subset_stringdb_ppi_df, file="../Data/subset_StringDB_interactions.csv",row.names=F, col.names=T, quote=F, sep=",")
