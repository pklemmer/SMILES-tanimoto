install.packages("RxnSim")#requires Rtools
install.packages("dplyr")
install.packages("utils")
packages <- c("readr","dplyr","utils", "RxnSim")
lapply(packages, library, character.only = TRUE)
Chemical_List_EPAPCS_2022_04_09 <- read_csv("C:/Chemical List EPAPCS-2022-04-09.csv")
df_smiles_unfiltered <- select(Chemical_List_EPAPCS_2022_04_09, "PREFERRED NAME", "SMILES")
df_smiles <- na.omit(df_smiles_unfiltered) #Cleaning up dataset by removing all entries without SMILES
df_smiles[, 'Similarity'] <- NA  #Adding a third column to df_smiles where the similarity values will be added upon calculation
View(df_smiles) #Why is a capital V required?
molA <- "[H][C@@]12COC3=C(C=C(OC)C(OC)=C3)[C@]1([H])C(=O)C1=CC=C3O[C@H](CC3=C1O2)C(C)=C" #This is the SMILES for Rotenone
for(row in 1:nrow(df_smiles)){calculation <- ms.compute(molA, df_smiles[row,"SMILES"], standardize = FALSE)}

