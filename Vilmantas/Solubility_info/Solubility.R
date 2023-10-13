install.packages("rvest")
library(rvest)
install.packages("httr")
library(httr)
install.packages("xml2")
library(xml2)
install.packages("jsonlite")
library(jsonlite)
install.packages("curl")
library(curl)
install.packages("RSelenium")
library(RSelenium)
install.packages("dplyr")
library(dplyr)

url <- "https://www.selleckchem.com/"
response <- GET(url)

# Example for extracting solubility information using CSS selectors
solubility_selector <- ".solubility > tbody > tr:nth-child(2) > td:nth-child(2)"
solubility <- html_text(html_nodes(response, solubility_selector))

# Define the compound names
compound_names <- c("Mitoxantrone", "Telaprevir", "Colchicine", "Benzydamine")

# Create an empty data frame to store the results
solubility_data <- data.frame(Compound = character(),
                              Solubility = character(),
                              stringsAsFactors = FALSE)

# Loop through each compound and search on the website
for (compound in compound_names) {
  # Construct the search URL for the compound
  encoded_compound <- URLencode(compound)
  search_url <- paste0("https://www.selleckchem.com/search.html?q=", encoded_compound)
  
  # Send a GET request to the search URL
  search_page <- read_html(search_url)
  
  # Extract the URL of the compound's page from the search results
  compound_url <- search_page %>%
    html_nodes(".medchem-product-name") %>%
    html_node("a") %>%
    html_attr("href")
  
  # Visit the compound's page
  compound_page <- read_html(compound_url)
  
  # Scrape the solubility information
  solubility <- compound_page %>%
    html_nodes(xpath = "//td[contains(text(), 'Solubility')]/following-sibling::td") %>%
    html_text() %>%
    trimws()
  
  # Add the compound and solubility to the data frame
  solubility_data <- rbind(solubility_data, data.frame(Compound = compound, Solubility = solubility))
}

# Print the solubility data
print(solubility_data)

--------------------
  library(httr)
library(rvest)

# Define the compound names
compound_names <- c("Mitoxantrone", "Telaprevir", "Colchicine", "Benzydamine")

# Create an empty data frame to store the results
solubility_data <- data.frame(Compound = character(),
                              Solubility = character(),
                              stringsAsFactors = FALSE)

# Loop through each compound and scrape solubility information
for (compound_name in compound_names) {
  # Form the compound's URL
  compound_url <- paste0("https://www.selleckchem.com/products/", URLencode(compound_name), ".html")
  
  # Send a GET request using httr
  response <- httr::GET(url = compound_url)
  
  # Parse the response content using rvest
  compound_page <- read_html(httr::content(response, as = "text"))
  
  # Extract solubility information using regex
  solubility_text <- compound_page %>%
    html_text() %>%
    gsub("\\s+", " ", .)  # Replace multiple spaces with single space
  
  solubility <- gsub("Solubility \\(25°C\\) In vitro Batch:", "", solubility_text)
  
  # Add the compound and solubility to the data frame
  solubility_data <- rbind(solubility_data, data.frame(Compound = compound_name, Solubility = solubility))
}

# Print the solubility data
print(solubility_data)








-------------------
  
  


# for pubchem
# Define the compound names
compound_names <- c(
  "Mitoxantrone", "Telaprevir", "Colchicine", "Benzydamine",
  "Varenicline Tartrate", "Naltrexone hydrochloride", "Rifapentine", "Capmatinib xHCl",
  "Tozasertib", "Iopromide", "Belumosudil", "Maraviroc", "Sotagliflozin", "Zoledronic Acid",
  "Decitabine", "Clopidogrel", "Abemaciclib", "Erdafitinib", "Ixazomib", "Methazolamide",
  "Dexmedetomidine", "Bepotastine Besilate", "Bromperidol", "Pirarubicin", "Nintedanib",
  "Cabazitaxel", "Raltegravir", "Sonidegib", "Eplerenone", "Ribociclib", "Carfilzomib",
  "Rolapitant hydrochloride", "D-Cycloserine", "Berbamine dihydrochloride", "Daphnetin",
  "Anethole trithione", "Avapritinib", "Galanthamine hydrobromide", "Racanisodamine",
  "Bilastine", "Selpercatinib", "Betulin", "Isavuconazole", "Sanguinarine chloride",
  "Spiperone", "Benfotiamine", "Raltitrexed", "DL-Adrenaline Hydrochloride", "Vismodegib",
  "Eltrombopag Olamine", "Doravirine", "Metformin", "Afatinib", "Reserpine", "Bedaquiline",
  "Sorafenib", "Imrecoxib", "Encorafenib", "Lusutrombopag", "Bictegravir", "Cobimetinib",
  "Topotecan hydrochloride", "Balsalazide sodium hydrate", "Grazoprevir", "Sulfaphenazole",
  "Sacubitril", "Desoximetasone", "Teneligliptin hydrobromide", "Tideglusib", "Pamabrom",
  "Pentamidine isethionate", "Gefitinib", "KPT330", "Rimegepant", "Indacaterol maleate",
  "Calcitriol", "Simeprevir", "Norgestrel", "Berberine", "Sotorasib"
)

# Create an empty data frame to store the results
solubility_data <- data.frame(Compound = character(),
                              Solubility = character(),
                              stringsAsFactors = FALSE)

# Loop through each compound and search on PubChem
for (compound in compound_names) {
  # Encode the compound name for the URL
  encoded_compound <- URLencode(compound)
  
  # Construct the PubChem search URL
  search_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", encoded_compound, "/JSON")
  
  # Send a GET request to the PubChem search URL
  search_response <- GET(search_url)
  
  # Parse the JSON response
  search_data <- fromJSON(content(search_response, "text"), flatten = TRUE)
  
  # Check if the compound was found
  if (length(search_data$PC_Compounds) > 0) {
    # Extract the CID (Compound Identifier) from the response
    cid <- search_data$PC_Compounds$CID[1]
    
    # Construct the PubChem compound summary URL
    summary_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/property/Solubility/JSON")
    
    # Send a GET request to the compound summary URL
    summary_response <- GET(summary_url)
    
    # Parse the JSON response
    summary_data <- fromJSON(content(summary_response, "text"), flatten = TRUE)
    
    # Check if solubility information is available
    if (!is.null(summary_data$PropertyTable$Properties$Solubility)) {
      # Extract the solubility value from the response
      solubility <- summary_data$PropertyTable$Properties$Solubility
    } else {
      # Set solubility to "Not found"
      solubility <- "Not found"
    }
    
    # Add the compound and solubility to the data frame
    solubility_data <- rbind(solubility_data, data.frame(Compound = compound, Solubility = solubility))
  } else {
    # Add a placeholder for compounds that were not found
    solubility_data <- rbind(solubility_data, data.frame(Compound = compound, Solubility = "Not found"))
  }
}

# Print the solubility data
print(solubility_data)
# WORKS! but not works

#From https://www.glpbio.com/
library(rvest)
library(httr)

# Define the compound names
compound_names <- c(
  "Mitoxantrone", "Telaprevir", "Colchicine", "Benzydamine",
  "Varenicline Tartrate", "Naltrexone hydrochloride", "Rifapentine", "Capmatinib xHCl",
  "Tozasertib", "Iopromide", "Belumosudil", "Maraviroc", "Sotagliflozin", "Zoledronic Acid",
  "Decitabine", "Clopidogrel", "Abemaciclib", "Erdafitinib", "Ixazomib", "Methazolamide",
  "Dexmedetomidine", "Bepotastine Besilate", "Bromperidol", "Pirarubicin", "Nintedanib",
  "Cabazitaxel", "Raltegravir", "Sonidegib", "Eplerenone", "Ribociclib", "Carfilzomib",
  "Rolapitant hydrochloride", "D-Cycloserine", "Berbamine dihydrochloride", "Daphnetin",
  "Anethole trithione", "Avapritinib", "Galanthamine hydrobromide", "Racanisodamine",
  "Bilastine", "Selpercatinib", "Betulin", "Isavuconazole", "Sanguinarine chloride",
  "Spiperone", "Benfotiamine", "Raltitrexed", "DL-Adrenaline Hydrochloride", "Vismodegib",
  "Eltrombopag Olamine", "Doravirine", "Metformin", "Afatinib", "Reserpine", "Bedaquiline",
  "Sorafenib", "Imrecoxib", "Encorafenib", "Lusutrombopag", "Bictegravir", "Cobimetinib",
  "Topotecan hydrochloride", "Balsalazide sodium hydrate", "Grazoprevir", "Sulfaphenazole",
  "Sacubitril", "Desoximetasone", "Teneligliptin hydrobromide", "Tideglusib", "Pamabrom",
  "Pentamidine isethionate", "Gefitinib", "KPT330", "Rimegepant", "Indacaterol maleate",
  "Calcitriol", "Simeprevir", "Norgestrel", "Berberine", "Sotorasib"
)

# Create an empty data frame to store the results
solubility_data <- data.frame(Compound = character(),
                              Solubility = character(),
                              stringsAsFactors = FALSE)

# Loop through each compound
for (compound in compound_names) {
  # Encode the compound name
  encoded_compound <- URLencode(compound)
  
  # Create the search URL
  search_url <- paste0("https://www.glpbio.com/search?name=", encoded_compound)
  
  # Fetch the HTML content
  response <- GET(search_url)
  
  # Check if the request was successful
  if (status_code(response) == 200) {
    # Parse the HTML content
    search_html <- read_html(content(response, "text"))
    
    # Extract the solubility information
    solubility <- search_html %>% html_nodes(".solubility") %>% html_text() %>% trimws()
    
    # Check if solubility information is found
    if (length(solubility) > 0) {
      solubility_data <- rbind(solubility_data, data.frame(Compound = compound, Solubility = solubility))
    } else {
      solubility_data <- rbind(solubility_data, data.frame(Compound = compound, Solubility = "Not found"))
    }
  } else {
    solubility_data <- rbind(solubility_data, data.frame(Compound = compound, Solubility = "Error"))
  }
}

# Print the solubility data
print(solubility_data)
