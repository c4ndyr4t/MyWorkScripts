install.packages("ggplot2")  # If not already installed
library(ggplot2)
data <- read.csv("decarboxylation.csv")
data$date <- as.Date(data$date)
data$Ratio <- as.numeric(data$Ratio)
ggplot(data, aes(x = date, y = Ratio) +
  geom_point() +
  geom_jitter() +
  labs(x = "Date", y = "Ratio") +
  theme_bw()
  
  data <- read.csv("decarboxylation.csv")
  head(data)
  
  install.packages("readxl")  # If not already installed
  library(readxl)
  library(ggplot2)
  data <- read_excel("decarboxylation.xlsx", sheet = 1)  # Replace "your_file.xlsx" with the actual file name and specify the sheet number or name if necessary
  data$Date <- as.Date(data$Date)
  ggplot(data, aes(x = Date, y = Ratio)) +
    geom_point() +
    geom_jitter() +
    labs(x = paste("Date", "(2022)"), y = paste("Ratio", "250/294")) +
    theme_bw()
  
  