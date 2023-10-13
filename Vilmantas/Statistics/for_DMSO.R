# Load the ggplot2 package
install.packages("ggplot2")
library(ggplot2)

# Create a data frame with your provided data
data <- data.frame(
  PlateID = c("PHD019086", "PHD019087", "PHD019088", "PHD019089", "PHD019090",
              "PHD019091", "PHD019092", "PHD019093", "PHD019094", "PHD019095",
              "PHD019096", "PHD019097", "PHD019098", "PHD019099", "PHD019100",
              "PHD019101", "PHD019102", "PHD019103"),
  FailureRate = c(0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7)
)

# Load the ggplot2 package
library(ggplot2)

# Create a data frame with your provided data
data <- data.frame(
  PlateID = c("PHD019086", "PHD019087", "PHD019088", "PHD019089", "PHD019090",
              "PHD019091", "PHD019092", "PHD019093", "PHD019094", "PHD019095",
              "PHD019096", "PHD019097", "PHD019098", "PHD019099", "PHD019100",
              "PHD019101", "PHD019102", "PHD019103"),
  FailureRate = c(0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7)
)

# Convert Plate IDs to numeric format
data$NumericPlateID <- as.numeric(substr(data$PlateID, 7, nchar(data$PlateID)))

# Create a bar chart with custom styling
ggplot(data, aes(x = reorder(PlateID, -FailureRate), y = FailureRate, fill = factor(FailureRate > 0))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("white", "#FF9999")) +  # Light red color
  labs(x = "Plate ID", y = "Failure Rate") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "#F4F4F4"),  # Light gray background
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)  # Rotate x-axis labels
  )

# Create a vector of failure rates
failure_rates <- c(0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7)

# Create a vector of plate IDs
plate_ids <- c("019086", "019087", "019088", "019089", "019090", "019091",
               "019092", "019093", "019094", "019095", "019096", "019097",
               "019098", "019099", "019100", "019101", "019102", "019103")

# Create a bar chart with red bars
barplot(failure_rates, names.arg = plate_ids, col = "red",
        xlab = "Failure Rate", ylab = "Plate ID", main = "Failure Rate by Plate ID")

# Create a vector of failure rates
failure_rates <- c(0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7)

# Create a vector of plate IDs
plate_ids <- c("019086", "019087", "019088", "019089", "019090", "019091",
               "019092", "019093", "019094", "019095", "019096", "019097",
               "019098", "019099", "019100", "019101", "019102", "019103")

# Create a bar chart with red bars and adjusted axis labels size
barplot(failure_rates, names.arg = plate_ids, col = "red",
        xlab = "Failure Rate", ylab = "Plate ID", main = "Failure Rate by Plate ID",
        cex.axis = 1)  # Adjust the size of axis labels

# Create a vector of failure rates
failure_rates <- c(0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7)




# Create a vector of failure rates
failure_rates <- c(0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7)

# Create a vector of plate IDs
plate_ids <- c("019086", "019087", "019088", "019089", "019090", "019091",
               "019092", "019093", "019094", "019095", "019096", "019097",
               "019098", "019099", "019100", "019101", "019102", "019103")

# Set up the plot with extra margin at the bottom
par(mar = c(7, 4, 4, 5) + 0.1)  # Increase the bottom margin

# Create a bar chart with red bars, rotated labels on the horizontal axis, and adjusted size
barplot(failure_rates, names.arg = plate_ids, col = "red",
        xlab = "Failure Rate", ylab = "Plate ID", main = "Failure Rate by Plate ID",
        las = 2, ylim = c(-1, max(failure_rates) + 1), cex.names = 0.8)

# Create a vector of failure rates
failure_rates <- c(0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7)

# Create a vector of plate IDs
plate_ids <- c("019086", "019087", "019088", "019089", "019090", "019091",
               "019092", "019093", "019094", "019095", "019096", "019097",
               "019098", "019099", "019100", "019101", "019102")

# Set up the plot with extra margin at the bottom
par(mar = c(7, 4, 4, 5) + 0.1)  # Increase the bottom margin

# Create a bar chart with custom colors, rotated labels, and improved aesthetics
barplot(failure_rates, names.arg = plate_ids, col = "skyblue",
        xlab = "Plate ID", ylab = "Failure Rate (PHD...)", main = "Compound (dissolved in DMSO) Failure Rate by Plate ID",
        las = 2, ylim = c(-1, max(failure_rates) + 1), cex.names = 0.7, border = "black")

# Add a horizontal grid for better readability
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")

# Add a legend
legend("topright", legend = c("Failure Rate"), fill = c("skyblue"))



#Another graph
# Load necessary libraries
install.packages("viridis")

library(ggplot2)
library(dplyr)

# Create a dataframe with your provided data (use NA for missing values)
data <- data.frame(
  Compound = c("Diphenylpyraline hydrochloride", "Eplerenone", "Pazopanib", 
               "Trichlormethiazide", "Cefdinir", "Fenoprofen calcium dihydrate", 
               "Thalidomide", "Betahistine mesylate", "Spiperone", 
               "Oxybuprocaine hydrochloride", "Gluconate Calcium", "Econazole", 
               "Zinc pyrithione", "Enoxacin"),
  Empirical_Formula = c("C19H23NO·HCl", "C24H30O6", "C21H23N7O2S", "C8H8Cl3N3O4S2", 
                        "C14H13N5O5S2", "C30H30CaO8", "C13H10N2O4", 
                        "C8H12N2·2(CH4O3S)", "C23H26FN3O2", "C17H29ClN2O3", 
                        "C12H22CaO14", "C18H15Cl3N2O", "C10H8N2O2S2Zn", "C15H17FN4O3"),
  Monoisotopic_Compound_Mass = c(282.35, 414.49, 437.52, 380.66, 395.03581, 
                                 558.63, 258.23, 328.4, 395.47, 308.20999, 
                                 430.37, 381.68, 317.7, 320.32),
  Experimental_Mass = c(79.0212, 79.0207, 79.0211, 79.0209, 79.0209, 
                        NA, 79.0208, NA, 79.021, 79.0209, NA, NA, 79.021, 79.021),
  DMSO_adduct = c(360.19973, 493.2256, 516.1849, NA, 474.0566, NA, 337.0847, NA, NA, NA, NA, NA, NA, 399.1498),
  DMSO_H_mass_accuracy_ppm = c(-0.126548303, -6.453963436, -1.392031329, 
                               -3.922997383, -3.922997383, NA, -5.188480409, NA, -2.657514356, 
                               -3.922997383, NA, NA, -2.657514356, -2.657514356),
  Intensity_DMSO_H = c(862972480, 545553620, 1063487552, 1138857728, 1150645120, 
                       NA, 1311306368, NA, 591469440, 1045218, NA, NA, 735952512, 1495830784),
  M_DMSO_H = c(4482251, 654922112, 68086831, NA, 899395, NA, 1126996, NA, NA, NA, NA, NA, NA, 224310656),
  M_H = c(2732518400, 308285280, 532161184, NA, 176196752, NA, 21898766, NA, 342052800, 
          1900751120, NA, NA, NA, 722989440),
  Comment = c("", "", "", "NH4", "", "neg ok", "MS2, sprst?", "both failed", "", "", "both failed", "both failed", "MS2, .molprob", ""),
  ratio_MH_C2H7SO = c(0.315815798, 1.769638888, 1.998431272, NA, 6.530455908, NA, 59.88037719, NA, 1.729175847, 0.000549897, NA, NA, NA, 2.068952465),
  ratio_MH_MC2H6SO = c(0.001640337, 2.1244028, 0.127944001, NA, 0.005104493, NA, 0.051463904, NA, NA, NA, NA, NA, NA, 0.310254402)
)

# Create the scatter plot
scatter_plot <- ggplot(data, aes(x = DMSO_H_mass_accuracy_ppm, 
                                 y = ratio_MH_C2H7SO, 
                                 color = ratio_MH_C2H7SO)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "green", high = "red") +
  labs(title = "DMSO Influence Analysis",
       x = "DMSO+H Mass Accuracy (ppm)",
       y = "Ratio [M+H]+ : [C2H7SO]+") +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

# Display the scatter plot
print(scatter_plot)


#ANOTHER ONE djKaled on the run

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create the bar chart
bar_chart <- ggplot(data, aes(x = reorder(Compound, Intensity_DMSO_H), 
                              y = Intensity_DMSO_H, 
                              fill = Compound)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Intensity of DMSO+H for Each Compound",
       x = "Compound",
       y = "Intensity of DMSO+H") +
  theme_minimal() +
  theme(legend.position="none")

# Display the bar chart
print(bar_chart)


