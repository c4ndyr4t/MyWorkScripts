library(stringr)
library(data.table)

# Function to recursively search for files with a specific extension
list_files <- function(path, extension) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  files <- files[str_detect(files, paste0("\\.", extension, "$"))]
  return(files)
}

# Function to extract keyword counts from a file
extract_keyword_counts <- function(file_path) {
  file_content <- fread(file_path)
  keyword_counts <- c(
    target_ions_not_found = sum(str_count(tolower(file_content), "target ions not found")),
    successful = sum(str_count(tolower(file_content), "successful"))
  )
  return(keyword_counts)
}

# Specify the root folder path
root_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/DATA/mzCloud/Data/Raw_data/Reruns"

# Get the paths of all experiment_status.txt files
file_paths <- list_files(root_folder, "txt")

# Initialize counters
target_ions_not_found_count <- 0
successful_count <- 0

# Iterate over each file and update the counters
for (file_path in file_paths) {
  counts <- extract_keyword_counts(file_path)
  target_ions_not_found_count <- target_ions_not_found_count + counts[["target_ions_not_found"]]
  successful_count <- successful_count + counts[["successful"]]
}

# Print the results
cat("Total occurrences of 'Target Ions not Found':", target_ions_not_found_count, "\n")
cat("Total occurrences of 'Successful':", successful_count, "\n")



root_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/DATA/mzCloud/Data/Raw_data/plate 12/Deimos/20220819_AnalytiCon_29_Rerun_1"
