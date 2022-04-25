  #Comments refer to the code above them.
  #This script generates Tanimoto/Jaccard chemical similarity scores. Here, it is used to determine what molecules in the dataset are most chemically similar to a reference chemical.
sessionInfo()
install.packages("RxnSim")
  #requires Rtools
install.packages("dplyr")
packages <- c("readr","dplyr","utils", "RxnSim", "ggplot2")
lapply(packages, library, character.only = TRUE)
  #The previous 2 lines make it easier to load packages and libraries without needing to enter library() every time

epapcs <- read_csv("Chemical List EPAPCS-2022-04-11.csv")
pc_eup <- read_csv("PubChem_compound_list_ElW0fZMu9pLBvH6l_N03iWB0JxQ494ui8YeQ7uqWgu_qj74.csv")
pc_agr <- read_csv("PubChem_compound_list_o-QFwS1hSN1_90ruyJYDwlQ_i1-vbTmhQ4Qi7ViVMOxYjAw.csv")
  #Importing datasets
epapcs <- select(epapcs, "PREFERRED NAME", "SMILES")
pc_eup <- select(pc_eup, "cmpdname", "isosmiles")
pc_agr <- select(pc_agr, "cmpdname", "isosmiles")
  #Selecting relevant columns

epapcs <- na.omit(epapcs)
pc_eup <- na.omit(pc_eup)
pc_agr <- na.omit(pc_agr)
  #Cleaning up dataset by removing all entries without SMILES

names(pc_eup)[1] <- "PREFERRED NAME"
names(pc_eup)[2] <- "SMILES"
names(pc_agr)[1] <- "PREFERRED NAME"
names(pc_agr)[2] <- "SMILES"
  #Renaming columns to be consistent across datasets
df_master <- rbind(epapcs, pc_eup, pc_agr)
df_master <- distinct(df_master)
  #Removing duplicates

molA <- "[H][C@@]12COC3=C(C=C(OC)C(OC)=C3)[C@]1([H])C(=O)C1=CC=C3O[C@H](CC3=C1O2)C(C)=C" 
  #This is the SMILES for Rotenone; defining molA as another SMILES value allows for comparisons of all kinds of pesticides.
df_master <- data.frame(lapply(df_master, as.character), stringsAsFactors = FALSE) 
  #For some reason, the df is seen as tables, not characters. This function forces all values to be characters
Tanimoto_coefficient <- vector("numeric",nrow(df_master))
  #Preparing a container for the calculated coefficient
for(row in 1:nrow(df_master)){
  calculation<- ms.compute(molA, df_master[row,"SMILES"], standardize = FALSE) 
  #This line repeats ms.compute for each row of df_smiles so that a comparison between Rotenone (=molA) and all other (available) pesticides (="SMILES") is made
  #RxnSim also allows comparisons between each pesticide in a list by using the ms.compute.sim.matrix function.
  Tanimoto_coefficient[row] <- calculation
}
df_master <- cbind (df_master,Tanimoto_coefficient)     
  #Tanimoto_coefficient is added to the df_smiles for easier viewing
df_master_sorted <- df_master[order(-df_master$Tanimoto_coefficient),]
View(df_master_sorted)
#write.csv(df_master_sorted, file = "Tanimoto-coefficient-Rotenone.csv")
  #Exports df_smiles_sorted to a .csv file with observations sorted in descending order according to their Tanimoto coefficient.
preferred_name <- df_master[,1]
p1 <- ggplot(data=df_master_sorted)+geom_col(mapping=aes(x= reorder(preferred_name,-Tanimoto_coefficient), y=Tanimoto_coefficient))
p1