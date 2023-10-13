install.packages("toast")
install.packages("robocopy")

#Copying all MS spectra

source_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/test_data_export/20220512_AnalytiCon_15_neg_2"
destination_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/Output_data_export"

# Function to recursively list all files with a specific extension
list_files_with_extension <- function(folder, extension) {
  files <- list.files(folder, full.names = TRUE, recursive = TRUE)
  file_list <- files[grepl(paste0("\\.", extension, "$"), files)]
  return(file_list)
}

# Check if a text file contains "MS1" or "MS2" (excluding "Ms1 -> Ms2")
check_spectra <- function(file) {
  text <- readLines(file)
  spectra_found <- any(grepl("MS1", text) & !grepl("Ms1 -> Ms2", text)) ||
    any(grepl("MS2", text) & !grepl("Ms1 -> Ms2", text))
  return(spectra_found)
}

# Check if a text file contains only MS1 spectra
check_only_ms1 <- function(file) {
  text <- readLines(file)
  ms1_found <- any(grepl("MS1", text) & !grepl("Ms1 -> Ms2", text))
  ms2_found <- any(grepl("MS2", text) & !grepl("Ms1 -> Ms2", text))
  return(ms1_found && !ms2_found)
}

# Get a list of all raw files in the source folder and its subfolders
raw_files <- list_files_with_extension(source_folder, "raw")
text_files <- list_files_with_extension(source_folder, "txt")

# Copy the raw files to the destination folder, excluding files with only MS1 spectra
for (raw_file in raw_files) {
  text_file <- gsub("\\.raw$", "_normal_log.txt", raw_file)
  if (text_file %in% text_files && check_spectra(text_file) && !check_only_ms1(text_file)) {
    destination_file <- file.path(destination_folder, basename(raw_file))
    if (!file.exists(destination_file)) {
      file.copy(raw_file, destination_folder)
    }
  }
}

#20230621

source_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/test_data_export/20220512_AnalytiCon_15_neg_2"
destination_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/Output_data_export"

# Function to recursively list all folders with a specific extension
list_folders_with_extension <- function(folder, extension) {
  folders <- list.dirs(folder, recursive = TRUE, full.names = TRUE)
  folder_list <- folders[grepl(paste0("\\.", extension, "$"), folders)]
  return(folder_list)
}

# Check if a folder contains raw files and meets the criteria for copying
check_folder <- function(folder) {
  raw_files <- list.files(folder, pattern = "\\.raw$", full.names = TRUE)
  text_files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
  ms1_found <- any(grepl("MS1", text_files) & !grepl("Ms1 -> Ms2", text_files))
  ms2_found <- any(grepl("MS2", text_files) & !grepl("Ms1 -> Ms2", text_files))
  return(length(raw_files) > 0 && ms1_found && !ms2_found)
}

# Copy the folders to the destination folder if they meet the criteria
source_folders <- list_folders_with_extension(source_folder, "raw")
for (folder in source_folders) {
  if (check_folder(folder)) {
    destination_subfolder <- gsub(source_folder, destination_folder, folder)
    system2("robocopy", args = c(folder, destination_subfolder, "/E"))
  }
}

# Send a notification once the script has completed
if (Sys.info()["sysname"] == "Windows") {
  command <- paste0("powershell -Command \"New-BurntToastNotification -Text 'Folder copying script has finished executing.'\"")
  system2(command)
}







#Lame

# Send a notification once the script has completed
if (Sys.info()["sysname"] == "Windows") {
  command <- paste0("powershell -Command \"New-BurntToastNotification -Text 'File copying script has finished executing.'\"")
  system2(command)
}

# 20230620 for also cCOMx files
# Function to recursively list all files with a specific extension
list_files_with_extension <- function(folder, extension) {
  files <- list.files(folder, full.names = TRUE, recursive = TRUE)
  file_list <- files[grepl(paste0("\\.", extension, "$"), files)]
  return(file_list)
}

# Check if a text file contains "MS1" or "MS2" (excluding "Ms1 -> Ms2")
check_spectra <- function(file) {
  text <- readLines(file)
  spectra_found <- any(grepl("MS1", text) & !grepl("Ms1 -> Ms2", text)) ||
    any(grepl("MS2", text) & !grepl("Ms1 -> Ms2", text))
  return(spectra_found)
}

# Check if a text file contains only MS1 spectra
check_only_ms1 <- function(file) {
  text <- readLines(file)
  ms1_found <- any(grepl("MS1", text) & !grepl("Ms1 -> Ms2", text))
  ms2_found <- any(grepl("MS2", text) & !grepl("Ms1 -> Ms2", text))
  return(ms1_found && !ms2_found)
}

# Get a list of all raw files in the source folder and its subfolders
raw_files <- list_files_with_extension(source_folder, "raw")
text_files <- list_files_with_extension(source_folder, "txt")
ccomx_files <- list_files_with_extension(source_folder, "ccomx")

# Create the "RAW" folder if it doesn't exist
raw_folder <- file.path(destination_folder, "RAW")
if (!file.exists(raw_folder)) {
  dir.create(raw_folder)
}

# Create the "cCOMx" folder if it doesn't exist
ccomx_folder <- file.path(destination_folder, "cCOMx")
if (!file.exists(ccomx_folder)) {
  dir.create(ccomx_folder)
}

# Copy the raw files to the "RAW" folder, excluding files with only MS1 spectra
for (raw_file in raw_files) {
  text_file <- gsub("\\.raw$", "_normal_log.txt", raw_file)
  if (text_file %in% text_files && check_spectra(text_file) && !check_only_ms1(text_file)) {
    destination_file <- file.path(raw_folder, basename(raw_file))
    if (!file.exists(destination_file)) {
      file.copy(raw_file, raw_folder)
    }
  }
}

# Copy the cCOMx files (excluding _unmerged.ccomx and MS1 spectra) to the "cCOMx" folder
for (ccomx_file in ccomx_files) {
  if (!grepl("_unmerged\\.ccomx$", ccomx_file)) {
    text_file <- gsub("\\.ccomx$", "_normal_log.txt", ccomx_file)
    if (text_file %in% text_files && check_spectra(text_file) && !check_only_ms1(text_file)) {
      destination_file <- file.path(ccomx_folder, basename(ccomx_file))
      if (!file.exists(destination_file)) {
        file.copy(ccomx_file, ccomx_folder)
      }
    }
  }
}

# 20230620 to pop out a table of waht was coppied and what was not

#Function to recursively list all files with a specific extension
list_files_with_extension <- function(folder, extension) {
  files <- list.files(folder, full.names = TRUE, recursive = TRUE)
  file_list <- files[grepl(paste0("\\.", extension, "$"), files)]
  return(file_list)
}

# Check if a text file contains "MS1" or "MS2" (excluding "Ms1 -> Ms2")
check_spectra <- function(file) {
  text <- readLines(file)
  spectra_found <- any(grepl("MS1", text) & !grepl("Ms1 -> Ms2", text)) ||
    any(grepl("MS2", text) & !grepl("Ms1 -> Ms2", text))
  return(spectra_found)
}

# Check if a text file contains only MS1 spectra
check_only_ms1 <- function(file) {
  text <- readLines(file)
  ms1_found <- any(grepl("MS1", text) & !grepl("Ms1 -> Ms2", text))
  ms2_found <- any(grepl("MS2", text) & !grepl("Ms1 -> Ms2", text))
  return(ms1_found && !ms2_found)
}

# Get a list of all raw files in the source folder and its subfolders
raw_files <- list_files_with_extension(source_folder, "raw")
text_files <- list_files_with_extension(source_folder, "txt")
ccomx_files <- list_files_with_extension(source_folder, "ccomx")

# Create the "RAW" folder if it doesn't exist
raw_folder <- file.path(destination_folder, "RAW")
if (!file.exists(raw_folder)) {
  dir.create(raw_folder)
}

# Create the "cCOMx" folder if it doesn't exist
ccomx_folder <- file.path(destination_folder, "cCOMx")
if (!file.exists(ccomx_folder)) {
  dir.create(ccomx_folder)
}

copied_files <- character()
not_copied_files <- character()

# Copy the raw files to the "RAW" folder, excluding files with only MS1 spectra
for (raw_file in raw_files) {
  text_file <- gsub("\\.raw$", "_normal_log.txt", raw_file)
  if (text_file %in% text_files && check_spectra(text_file) && !check_only_ms1(text_file)) {
    destination_file <- file.path(raw_folder, basename(raw_file))
    if (!file.exists(destination_file)) {
      file.copy(raw_file, raw_folder)
      copied_files <- append(copied_files, raw_file)
    } else {
      not_copied_files <- append(not_copied_files, raw_file)
    }
  }
}

# Copy the cCOMx files (excluding _unmerged.ccomx and MS1 spectra) to the "cCOMx" folder
for (ccomx_file in ccomx_files) {
  if (!grepl("_unmerged\\.ccomx$", ccomx_file)) {
    text_file <- gsub("\\.ccomx$", "_normal_log.txt", ccomx_file)
    if (text_file %in% text_files && check_spectra(text_file) && !check_only_ms1(text_file)) {
      destination_file <- file.path(ccomx_folder, basename(ccomx_file))
      if (!file.exists(destination_file)) {
        file.copy(ccomx_file, ccomx_folder)
        copied_files <- append(copied_files, ccomx_file)
      } else {
        not_copied_files <- append(not_copied_files, ccomx_file)
      }
    }
  }
}

# Print the summary table
if (length(copied_files) > 0) {
  summary_table <- data.frame("Copied Files" = copied_files)
  print(summary_table)
} else {
  print("No files were copied.")
}

# Print the list of not copied files
if (length(not_copied_files) > 0) {
  print("Not Copied Files:")
  for (file in not_copied_files) {
    print(file)
  }
}

#20230620 whole folder
source_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/test_data_export/20220512_AnalytiCon_15_neg_2"
destination_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/Output_data_export"

# Function to recursively list all files with a specific extension
list_files_with_extension <- function(folder, extension) {
  files <- list.files(folder, full.names = TRUE, recursive = TRUE)
  file_list <- files[grepl(paste0("\\.", extension, "$"), files)]
  return(file_list)
}

# Check if a text file contains "MS1" or "MS2" (excluding "Ms1 -> Ms2")
check_spectra <- function(file) {
  text <- readLines(file)
  spectra_found <- any(grepl("MS1", text) & !grepl("Ms1 -> Ms2", text)) ||
    any(grepl("MS2", text) & !grepl("Ms1 -> Ms2", text))
  return(spectra_found)
}

# Check if a text file contains only MS1 spectra
check_only_ms1 <- function(file) {
  text <- readLines(file)
  ms1_found <- any(grepl("MS1", text) & !grepl("Ms1 -> Ms2", text))
  ms2_found <- any(grepl("MS2", text) & !grepl("Ms1 -> Ms2", text))
  return(ms1_found && !ms2_found)
}

# Get a list of all raw files in the source folder and its subfolders
raw_files <- list_files_with_extension(source_folder, "raw")
text_files <- list_files_with_extension(source_folder, "txt")

# Copy the folders containing the raw files to the destination folder
for (raw_file in raw_files) {
  folder <- dirname(raw_file)
  text_file <- gsub("\\.raw$", "_normal_log.txt", raw_file)
  if (text_file %in% text_files && check_spectra(text_file) && !check_only_ms1(text_file)) {
    destination_folder_path <- gsub(source_folder, destination_folder, folder)
    if (!dir.exists(destination_folder_path)) {
      dir.create(destination_folder_path, recursive = TRUE)
    }
    if (!file.exists(file.path(destination_folder_path, basename(raw_file)))) {
      file.copy(raw_file, destination_folder_path)
    }
  }
}

# Send a notification once the script has completed
if (Sys.info()["sysname"] == "Windows") {
  command <- paste0("powershell -Command \"New-BurntToastNotification -Text 'Folder copying script has finished executing.'\"")
  system2(command)
}

#20230621 one better raw file the other, - less good
source_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/test_data_export/20220512_AnalytiCon_15_neg_2"
destination_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/Output_data_export"

# Get a list of all raw files in the source folder and its subfolders
raw_files <- list.files(source_folder, pattern = "\\.raw$", full.names = TRUE, recursive = TRUE)

# Identify duplicate raw files based on file name
duplicate_files <- raw_files[duplicated(basename(raw_files))]

# Function to extract MS level from the "_normal_log.txt" file
get_ms_level <- function(file) {
  txt_file <- gsub("\\.raw$", "_normal_log.txt", file)
  if (file.exists(txt_file)) {
    text <- readLines(txt_file)
    ms_level <- sub(".*MS(\\d+).*", "\\1", text)
    return(ms_level)
  } else {
    return(NA)
  }
}

# Function to get file size in MB
get_file_size <- function(file) {
  file_info <- file.info(file)
  file_size <- file_info$size / (1024^2)  # Convert bytes to MB
  return(file_size)
}

# Copy the duplicate raw files to the destination folder while preserving the directory structure
for (file in duplicate_files) {
  destination_file <- file.path(destination_folder, file)
  if (!file.exists(destination_file)) {
    dir.create(dirname(destination_file), recursive = TRUE)
    file.copy(file, destination_file)
  }
}

# Store the results in a data frame
results <- data.frame(
  Raw_File = basename(duplicate_files),
  MS_Level = sapply(duplicate_files, get_ms_level),
  File_Size_MB = sapply(duplicate_files, get_file_size)
)

# Print the results
print(results)

#20230704Extracting information from a provided list (for dmso analysis)
source_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/DATA/mzCloud/Data/Raw_data/Target Mol"
destination_folder <- "//ltvil-freenas5.thermofisher.lt/PTVG_Data/Vilmantas/Output_data_export"

# List of IDs to search for
ids <- c(
  "T3098", "T2330", "T3489", "T6552", "T1011", "T2850", "T6218", "T1562", "T2923", "T0213",
  "T1109", "T1330", "T1639", "T2329", "T6214", "T6211", "T0717L", "T6708", "T2384", "T2385",
  "T1575", "T1465", "T1546", "T4093", "T6588", "T6190", "T6867", "T14830", "T6628", "T2543",
  "T0083", "T3724", "T1589", "T2851", "T5109", "T8222", "T0280", "T6825", "TQ0223", "T1181",
  "T4686", "T6115", "T3112", "T8541", "T11897", "T13693", "T17226", "T9251", "T19862", "T0378",
  "T2115", "T1797", "T6968", "T7903", "T1747"
)

# Function to search and copy files recursively
search_and_copy_files <- function(ids, source_folder, destination_folder) {
  # Initialize a list to store missing IDs
  missing_ids <- vector("character")
  
  # Recursive function to search for files matching the pattern
  search_files_recursive <- function(folder, file_pattern) {
    files <- list.files(folder, pattern = file_pattern, full.names = TRUE)
    subfolders <- list.dirs(folder, recursive = FALSE)
    
    matching_files <- files
    
    if (length(subfolders) > 0) {
      for (subfolder in subfolders) {
        matching_files <- c(matching_files, search_files_recursive(subfolder, file_pattern))
      }
    }
    
    return(matching_files)
  }
  
  # Iterate over the IDs
  for (id in ids) {
    # Generate the file name pattern to search for
    file_pattern <- paste0(id, "[-](pos|neg)\\.ccomx")
    
    # Search for files matching the pattern recursively
    matching_files <- search_files_recursive(source_folder, file_pattern)
    
    # Check if any files are found
    if (length(matching_files) > 0) {
      # Copy the files to the destination folder
      for (file in matching_files) {
        destination_file <- file.path(destination_folder, basename(file))
        if (!file.exists(destination_file)) {
          file.copy(file, destination_folder)
        }
      }
    } else {
      # Add the missing ID to the list
      missing_ids <- c(missing_ids, id)
    }
  }
  
  # Print the missing IDs
  if (length(missing_ids) > 0) {
    cat("Missing IDs:")
    cat(missing_ids, sep = ", ")
    cat("\n")
  }
}

# Perform the search and copy
search_and_copy_files(ids, source_folder, destination_folder)










