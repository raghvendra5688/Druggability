library(data.table)
library(ggplot2)

setwd("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/Druggability_Score/scripts/")

convert_to_listed_df <- function(df, column="", rename_column="")
{
  unique_drug_ids <- unique(df$drug_id)
  print(length(unique_drug_ids))
  revised_df <- NULL
  for (i in 1:length(unique_drug_ids))
  {
    drug_id <- unique_drug_ids[i]
    values <- df[df$drug_id==drug_id,column]
    if (length(values)>1)
    {
      revised_values <- paste0(values,collapse=" ; ")
    }else{
      revised_values <- values
    }
    revised_df <- rbind(revised_df,cbind(drug_id, revised_values))
  }
  revised_df <- as.data.frame(revised_df)
  colnames(revised_df) <- c("drug_id",rename_column)
  return(revised_df)
}

# #Load the drug database 
# ################################################################################
# load("Data/drug_data.RData")
# drug_general <- as.data.frame(DRUGBANK$general)
# colnames(drug_general)[1] <- "drug_id"
# drug_status <- as.data.frame(DRUGBANK$status)
# rev_drug_status <- convert_to_listed_df(drug_status, column="group", rename_column="collated_group")
# drug_df <- merge(drug_general,rev_drug_status, by="drug_id",all.x=T)
# drug_category <- as.data.frame(DRUGBANK$category)
# colnames(drug_category)[1] <- "drug_id"
# rev_drug_category <- convert_to_listed_df(drug_category,column="category", rename_column="collated_category")
# drug_df <- merge(drug_df,rev_drug_category, by="drug_id",all.x=T)
# 
# #Remove text related to mechanism of action, protein binding and toxicity and keep as separate dataframe
# ################################################################################
# drug_moa_df <- drug_df[,c("drug_id","name","description","indication","pharmacodynamics","mechanism_of_action","protein_binding","toxicity")]
# write.table(drug_moa_df, file="Data/drug_MOA_Description.csv",row.names=F, col.names=T, quote=F, sep=" | ")
# 
# drug_df <- drug_df[,!(colnames(drug_df) %in% c("description","indication","pharmacodynamics","mechanism_of_action","protein_binding","toxicity"))]

#Get the drugbank dataset with information with atleast one transporter, target or enzyme (human / non-human)
################################################################################
drug_basic_df <- fread("../Data/drugs.csv",header=T)
drug_basic_df <- as.data.frame(drug_basic_df)

#Get the drugbank dataset with calculated properties
################################################################################
drug_smiles_df <- fread("../Data/smiles.csv",header=T,sep="$")
drug_smiles_df <- as.data.frame(drug_smiles_df)

#Perform inner join and keep those drugs whose SMILES are available and they have one target (human / non-human)
############################################################################################
rev_drug_drugs_df <- merge(x = drug_smiles_df, y =  drug_basic_df, by = "drug_id")

#Get drugs interaction with human protein targets, transporters and enzymes
################################################################################
drug_target_df <- fread("../Data/drug2target_human.csv",header=T)
drug_target_df <- as.data.frame(drug_target_df)
drug_target_df$target_type <- "protein"
rev_drug_target_df <- drug_target_df[,c(1,2,3,6,4,7)]
colnames(rev_drug_target_df) <- c("drug_id","partner_id","gene_name","inducer","inhibitor","target_type")
drug_enzyme_df <- fread("../Data/drug2enzyme_human.csv",header=T)
drug_enzyme_df <- as.data.frame(drug_enzyme_df)
drug_enzyme_df$target_type <- "enzyme"
rev_drug_enzyme_df <- drug_enzyme_df[,c(1,2,3,5,6,7)]
drug_transporter_df <- fread("../Data/drug2transporter_human.csv",header=T)
drug_transporter_df <- as.data.frame(drug_transporter_df)
drug_transporter_df$target_type <- "transporter"
rev_drug_transporter_df <- drug_transporter_df[,c(1,2,3,5,6,7)]

drug_drug_protein_df <- as.data.frame(rbind(rev_drug_target_df, rev_drug_enzyme_df, rev_drug_transporter_df))
drug_drug_protein_df <- unique(drug_drug_protein_df)

#Make list version drug_id, partner_id
drug_drugid_partnerid_df <- convert_to_listed_df(drug_drug_protein_df, column="partner_id",rename_column = "collated_partner_id")
#Make list version drug_id, gene_name
drug_drugid_gene_names_df <- convert_to_listed_df(drug_drug_protein_df, column="gene_name",rename_column = "collated_gene_name")
#Make list version drug_id, target_type
drug_drugid_target_type_df <- convert_to_listed_df(drug_drug_protein_df, column="target_type",rename_column = "collated_target_type")
#Make list version drug_id, inhibitor
drug_drugid_inhibitor_df <- convert_to_listed_df(drug_drug_protein_df, column="inhibitor",rename_column = "collated_inhibitor")
#Make list version drug_id, inhibitor
drug_drugid_inducer_df <- convert_to_listed_df(drug_drug_protein_df, column="inducer",rename_column = "collated_inducer")

#Combine all drugbank protein interaction information
drug_drug_protein_interaction_df <- as.data.frame(cbind(drug_drugid_partnerid_df,drug_drugid_gene_names_df$collated_gene_name,drug_drugid_target_type_df$collated_target_type, drug_drugid_inhibitor_df$collated_inhibitor, drug_drugid_inducer_df$collated_inducer))
colnames(drug_drug_protein_interaction_df) <- c("drug_id","collated_partner_id","collated_gene_name","collated_target_type","collated_inhibitor","collated_inducer")

#Merge the drugbank information with drugbank protein interaction information
big_drug_data_df <- merge(x=rev_drug_drugs_df, y = drug_drug_protein_interaction_df)
print(dim(big_drug_data_df))
write.table(big_drug_data_df, file="../Data/drug_protein_interactions.csv",row.names=F, col.names=T, quote=F, sep="|")
