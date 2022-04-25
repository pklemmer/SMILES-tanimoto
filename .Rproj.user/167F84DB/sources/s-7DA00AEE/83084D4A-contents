library(readr)
pc_eup <- read_csv("PubChem_compound_list_ElW0fZMu9pLBvH6l_N03iWB0JxQ494ui8YeQ7uqWgu_qj74.csv")
pc_agr <- read_csv("PubChem_compound_list_o-QFwS1hSN1_90ruyJYDwlQ_i1-vbTmhQ4Qi7ViVMOxYjAw.csv")
install.packages("dplyr")
library(dplyr)

pc_eup <- select(pc_eup, "cmpdname", "isosmiles")
pc_agr <- select(pc_agr, "cmpdname", "isosmiles")
  #Selecting the name and SMILE column from the df
names(pc_eup)[1] <- "PREFERRED NAME"
names(pc_eup)[2] <- "SMILES"
names(pc_agr)[1] <- "PREFERRED NAME"
names(pc_agr)[2] <- "SMILES"
  #Renaming columns to be consistent across datasets
View(pc_eup)
View(pc_agr)
