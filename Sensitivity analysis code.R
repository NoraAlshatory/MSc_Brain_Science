
# code for the specificity and sensitivity between PRESS and STEAM at TE =35, 78 and 97ms.
# indicates which sequence is the best in theses TEs. 



# Load required packages

if (!require(pROC)) install.packages("pROC", dependencies=TRUE)
library(readr)
library(dplyr)
library(pROC)

# Load the data
press_data <- read_csv("press.csv")
steam_data <- read_csv("steam .csv")

# Clean and structure the STEAM data by removing duplicated metabolites
steam_data_cleaned <- steam_data %>%
  slice(2:n()) %>%
  select(TE, `0mM at 35`, `1mM at 35`, `2mM at 35`, `4mM at 35`, `8mM at 35`, 
         `0mM at 78`, `1mM at 78`, `2mM at 78`, `4mM at 78`, `8mM at 78`, 
         `0mM at 97`, `1mM at 97`, `2mM at 97`, `4mM at 97`, `8mM at 97`) %>%
  distinct()

# Rename columns for clarity
colnames(steam_data_cleaned) <- c("Metabolite", "STEAM_Conc_35", "STEAM_Conc_35_1", "STEAM_Conc_35_2", "STEAM_Conc_35_4", "STEAM_Conc_35_8",
                                  "STEAM_Conc_78", "STEAM_Conc_78_1", "STEAM_Conc_78_2", "STEAM_Conc_78_4", "STEAM_Conc_78_8",
                                  "STEAM_Conc_97", "STEAM_Conc_97_1", "STEAM_Conc_97_2", "STEAM_Conc_97_4", "STEAM_Conc_97_8")

# Extract relevant columns and rows for PRESS
press_data_cleaned <- press_data %>%
  slice(2:n()) %>%
  select(TE, `0mM at 35`, `1mM at 35`, `2mM at 35`, `4mM at 35`, `8mM at 35`, 
         `0mM at 78`, `1mM at 78`, `2mM at 78`, `4mM at 78`, `8mM at 78`, 
         `0mM at 97`, `1mM at 97`, `2mM at 97`, `4mM at 97`, `8mM at 97`) %>%
  distinct()

# Rename columns for clarity
colnames(press_data_cleaned) <- c("Metabolite", "PRESS_Conc_35", "PRESS_Conc_35_1", "PRESS_Conc_35_2", "PRESS_Conc_35_4", "PRESS_Conc_35_8",
                                  "PRESS_Conc_78", "PRESS_Conc_78_1", "PRESS_Conc_78_2", "PRESS_Conc_78_4", "PRESS_Conc_78_8",
                                  "PRESS_Conc_97", "PRESS_Conc_97_1", "PRESS_Conc_97_2", "PRESS_Conc_97_4", "PRESS_Conc_97_8")

# Convert columns to numeric
press_data_cleaned[, -1] <- lapply(press_data_cleaned[, -1], as.numeric)
steam_data_cleaned[, -1] <- lapply(steam_data_cleaned[, -1], as.numeric)

# Check for NA values and handle them if necessary
press_data_cleaned <- na.omit(press_data_cleaned)
steam_data_cleaned <- na.omit(steam_data_cleaned)

# Ensure both dataframes have the same number of rows
stopifnot(nrow(press_data_cleaned) == nrow(steam_data_cleaned))

# Create condition labels (0 for control, 1 for disease)
# Assuming equal control and disease samples for illustration
condition <- rep(c(0, 1), length.out = nrow(press_data_cleaned))

# Combine condition labels with the cleaned data
press_data_cleaned$Condition <- condition
steam_data_cleaned$Condition <- condition

# Perform ROC analysis for each TE and sequence
results <- data.frame()
TEs <- c("35", "78", "97")

for (TE in TEs) {
  for (metabolite in c("2HG", "Cysta")) {
    # Perform ROC analysis for PRESS
    suppressMessages({
      suppressWarnings({
        roc_press <- roc(press_data_cleaned$Condition, press_data_cleaned[[paste0("PRESS_Conc_", TE)]])
      })
    })
    auc_press <- auc(roc_press)
    ci_press <- ci.auc(roc_press)
    
    # Calculate sensitivity and specificity for the best threshold
    coords_press <- coords(roc_press, "best", ret = c("sensitivity", "specificity"))
    sensitivity_press <- coords_press["sensitivity"]
    specificity_press <- coords_press["specificity"]
    
    # Add results to the summary table for PRESS
    results <- rbind(results, data.frame(
      TE = TE,
      Method = "PRESS",
      Biomarker = metabolite,
      Sensitivity = sensitivity_press * 100,
      Specificity = specificity_press * 100,
      AUC = auc_press,
      `95% CI` = paste0(round(ci_press[1], 4), "-", round(ci_press[3], 4))
    ))
    
    # Perform ROC analysis for STEAM
    suppressMessages({
      suppressWarnings({
        roc_steam <- roc(steam_data_cleaned$Condition, steam_data_cleaned[[paste0("STEAM_Conc_", TE)]])
      })
    })
    auc_steam <- auc(roc_steam)
    ci_steam <- ci.auc(roc_steam)
    
    # Calculate sensitivity and specificity for the best threshold
    coords_steam <- coords(roc_steam, "best", ret = c("sensitivity", "specificity"))
    sensitivity_steam <- coords_steam["sensitivity"]
    specificity_steam <- coords_steam["specificity"]
    
    # Add results to the summary table for STEAM
    results <- rbind(results, data.frame(
      TE = TE,
      Method = "STEAM",
      Biomarker = metabolite,
      Sensitivity = sensitivity_steam * 100,
      Specificity = specificity_steam * 100,
      AUC = auc_steam,
      `95% CI` = paste0(round(ci_steam[1], 4), "-", round(ci_steam[3], 4))
    ))
  }
}

# Print the results
print(results)

# Plot ROC curves for each TE, method, and biomarker
par(mfrow = c(length(TEs), 1)) # Ensures there are plots in a column

for (TE in TEs) {
  for (metabolite in c("2HG", "Cysta")) {
    suppressMessages({
      suppressWarnings({
        roc_press <- roc(press_data_cleaned$Condition, press_data_cleaned[[paste0("PRESS_Conc_", TE)]])
        roc_steam <- roc(steam_data_cleaned$Condition, steam_data_cleaned[[paste0("STEAM_Conc_", TE)]])
      })
    })
    
    plot(roc_press, col = "darkorange", main = paste("ROC Curves for TE =", TE, "and Metabolite =", metabolite), lwd = 2)
    plot(roc_steam, col = "blue", add = TRUE, lwd = 2)
    legend("bottomright", legend = c(paste("PRESS (AUC =", round(auc(roc_press), 2), ")"),
                                     paste("STEAM (AUC =", round(auc(roc_steam), 2), ")")),
           col = c("darkorange", "blue"), lwd = 2)
  }
}



