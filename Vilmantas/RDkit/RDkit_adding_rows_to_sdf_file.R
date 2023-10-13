install.packages("ChemmineR")
install.packages("rdkit")
library(ChemmineR)

install.packages("remotes")
library(remotes)
install_version("ChemmineR", version = "3.4.0")

install.packages("rcdk")
library(rcdk)
sdf_data <- parse.sdf("L1000_compound_library_revVP.sdf")
