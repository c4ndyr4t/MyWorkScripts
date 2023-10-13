install.packages("httr")
install.packages("xml2")
library(httr)
library(xml2)

getCompoundInfo <- function(compound) {
  # Prepare the API URL
  base_url <- "https://go.drugbank.com"
  endpoint <- "/structures/smiles"
  query <- list(smiles = compound)
  
  # Send a GET request to the DrugBank API
  response <- GET(url = paste0(base_url, endpoint), query = query)
  
  # Check if the request was successful
  if (status_code(response) == 200) {
    # Extract the XML content from the response
    content_xml <- content(response, as = "text")
    
    # Parse the XML data
    xml_data <- read_xml(content_xml)
    
    return(xml_data)
  } else {
    # Return NULL if there was an error
    return(NULL)
  }
}

getCompoundGroup <- function(xml_data) {
  # Find the relevant XML nodes containing group information
  group_nodes <- xml_find_all(xml_data, "//groups/group")
  
  # Extract the group categories
  group_categories <- xml_attr(group_nodes, "category")
  
  # Check if the compound is approved, vet approved, or illicit
  if ("approved" %in% group_categories) {
    return("Approved")
  } else if ("vet_approved" %in% group_categories) {
    return("Vet Approved")
  } else if ("illicit" %in% group_categories) {
    return("Illicit")
  } else {
    return("Unknown")
  }
}

compound <- "Sincalide"

# Get compound information from DrugBank
compound_info <- getCompoundInfo(compound)

if (!is.null(compound_info)) {
  # Get the compound's group category
  group_category <- getCompoundGroup(compound_info)
  
  cat("Compound:", compound, "\n")
  cat("Group Category:", group_category, "\n")
} else {
  cat("Failed to retrieve compound information from DrugBank.\n")
}

# 20230622 PubChem 

install.packages("rentrez")
install.packages("jsonlite")
library(rentrez)
library(jsonlite)

getCompoundClassification <- function(compound) {
  # Search for the compound in PubChem
  search_result <- entrez_search(db = "pccompound", term = compound)
  
  if (search_result$count > 0) {
    # Get the PubChem ID of the first result
    pubchem_id <- search_result$ids[1]
    
    # Fetch the compound record from PubChem
    compound_record <- entrez_fetch(db = "pccompound", id = pubchem_id, rettype = "json")
    
    # Parse the JSON data
    compound_data <- fromJSON(compound_record, flatten = TRUE)
    
    # Extract the classification information
    classification <- compound_data$PC_Compounds$PC_Compound$props$IUPACName
    # You can extract other relevant information as needed
    
    return(classification)
  } else {
    # Return NULL if the compound is not found
    return(NULL)
  }
}

compound <- "CID1983"

# Get compound classification from PubChem
classification <- getCompoundClassification(compound)

if (!is.null(classification)) {
  cat("Compound:", compound, "\n")
  cat("Classification:", classification, "\n")
} else {
  cat("Compound not found in PubChem.\n")
}


