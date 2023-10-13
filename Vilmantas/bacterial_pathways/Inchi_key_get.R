install.packages("pubchem")
library(pubchem)
install.packages("rentrez")
library(rentrez)
install.packages("xml2")
library(xml2)
install.packages("rcdk")
library(rcdk)
# Install required packages
install.packages(c("readxl", "httr"))

# Load required packages
library(readxl)
library(httr)

install.packages("openxlsx")
library(openxlsx)

install.packages("writexl")
library(writexl)

# Load the required library
install.packages(xlsx)
library(xlsx)

#Remove environment
rm(list = ls())




#INCHI GETTER for one compound

# Define the compound name or identifier
compound_name <- "1-butanol"


# Make a query to the PubChem REST API using compound name
response <- httr::GET(
  url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", compound_name, "/property/InChI/JSON")
)

# Check if the request was successful
if (httr::status_code(response) == 200) {
  # Extract the InChI string from the response
  data <- httr::content(response)
  inchi_string <- data$PropertyTable$Properties[[1]]$InChI
  
  # Print the InChI string
  cat("InChI String:", inchi_string, "\n")
} else {
  cat("Error:", httr::status_code(response), "\n")
}

#INCHI, string and key and SMILES getter for mulittutde of compounds

# Read chemical names from Excel file
excel_file <- "C:/Users/vilmantas.pedisius/OneDrive - Thermo Fisher Scientific/Desktop/chemical_names.xlsx"  # Replace with your file path
sheet_name <- "MG1655_list"  # Replace with the sheet name containing the names
chemical_names <- readxl::read_excel(excel_file, sheet = sheet_name, col_names = FALSE)[[1]]
  
# Function to fetch InChI, InChI Key, and SMILES for a given chemical name
fetch_chemical_info <- function(chemical_name) {
  encoded_name <- URLencode(chemical_name)
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", encoded_name, "/property/InChI,InChIKey,CanonicalSMILES/JSON")
  response <- httr::GET(url)
  
  if (httr::status_code(response) == 200) {
    data <- httr::content(response)
    if ("Properties" %in% names(data$PropertyTable)) {
      inchi <- data$PropertyTable$Properties[[1]]$InChI
      inchikey <- data$PropertyTable$Properties[[1]]$InChIKey
      smiles <- data$PropertyTable$Properties[[1]]$CanonicalSMILES
      return(c(inchi, inchikey, smiles))
    } else {
      return(c(NA, NA, NA))
    }
  } else {
    return(c(NA, NA, NA))
  }
}

# Apply the function to fetch chemical information for each name
chemical_info <- lapply(chemical_names, fetch_chemical_info)

# Convert the list to a data frame
chemical_info_df <- data.frame(ChemicalName = chemical_names,
                               InChI = sapply(chemical_info, "[[", 1),
                               InChIKey = sapply(chemical_info, "[[", 2),
                               SMILES = sapply(chemical_info, "[[", 3),
                               stringsAsFactors = FALSE)

# Save the chemical information to a separate Excel file
output_file <- "C:/Users/vilmantas.pedisius/OneDrive - Thermo Fisher Scientific/Desktop/chemical_info.xlsx"
write_xlsx(chemical_info_df, path = output_file, col_names = TRUE)

#Inchi Keys to Names Inchies strings and smiles

install.packages("future.apply")
install.packages("httr")
# Load the required libraries
library(httr)
library(future.apply)

install.packages("readxl")
library(readxl)

# Set up parallel processing
plan(multisession)

# Read chemical names from Excel file
excel_file <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/Smilestoinchikey/chemical_info_Vytas_.xlsx"  # Replace with your file path
sheet_name <- "All_KEYS"  # Replace with the sheet name containing the names
chemical_names <- readxl::read_excel(excel_file, sheet = sheet_name, col_names = FALSE)[[1]]

# Function to fetch chemical information based on InChI Key
fetch_chemical_info <- function(inchikey) {
  encoded_key <- URLencode(inchikey)
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", encoded_key, "/property/IUPACName,InChI,CanonicalSMILES/JSON")
  response <- httr::GET(url)
  
  if (httr::status_code(response) == 200) {
    data <- httr::content(response)
    if ("Properties" %in% names(data$PropertyTable)) {
      name <- data$PropertyTable$Properties[[1]]$IUPACName
      inchi <- data$PropertyTable$Properties[[1]]$InChI
      smiles <- data$PropertyTable$Properties[[1]]$CanonicalSMILES
      return(c(name, inchi, smiles))
    } else {
      return(c(NA, NA, NA))
    }
  } else {
    return(c(NA, NA, NA))
  }
}

# Apply the function to fetch chemical information for each InChI Key
chemical_info <- lapply(chemical_names, fetch_chemical_info)

# Convert the list to a data frame
chemical_info_df <- data.frame(InchiKey = chemical_names,
                               Name = sapply(chemical_info, function(x) if (length(x) >= 1) x[[1]] else NA),
                               InChIString = sapply(chemical_info, function(x) if (length(x) >= 2) x[[2]] else NA),
                               SMILES = sapply(chemical_info, function(x) if (length(x) >= 3) x[[3]] else NA),
                               stringsAsFactors = FALSE)

# Save the chemical information to a separate Excel file
output_file <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/Smilestoinchikey/chemical_info_output.xlsx"
openxlsx::write.xlsx(chemical_info_df, file = output_file, sheetName = "ChemicalInfo", rowNames = FALSE)
