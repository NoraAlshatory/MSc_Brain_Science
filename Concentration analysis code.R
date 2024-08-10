

# Concentration analysis for 6ms at STEAM at 8mM.
#ATTH2O and WCONC






# Load necessary library
library(ggplot2)
library(gridExtra)
library(grid)  # This library provides the textGrob function

# Data for 2HG
data_2HG <- data.frame(
  Phantom = c(1, 2, 3, 4, 5),
  Observed_Concentration = c(0.156, 0, 0.32, 10.756, 16.608),
  Truth_Concentration = c(0, 1, 2, 4, 8)
)

# Perform the linear regression for 2HG
regression_model_2HG <- lm(Observed_Concentration ~ Truth_Concentration, data = data_2HG)
regression_summary_2HG <- summary(regression_model_2HG)

# Plot for 2HG
plot_2HG <- ggplot(data_2HG, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
  labs(title = "Actual vs Observed Concentrations of 2HG",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 8) +
  ylim(0, 20) +
  annotate("text", x = 3, y = 18, label = paste("R² =", round(regression_summary_2HG$r.squared, 2), 
                                                "\nIntercept =", round(regression_model_2HG$coefficients[1], 2), 
                                                "\np-value =", round(regression_summary_2HG$coefficients[2, 4], 2)),
           size = 5, hjust = 0, vjust = 1, color = "blue") +
  theme(panel.border = element_rect(color = "gray", fill = NA, size = 1),
        plot.background = element_rect(color = "gray", fill = NA, size = 1),
        axis.line = element_line(color = "gray"))

# Data for Cysta
data_Cysta <- data.frame(
  Phantom = c(1, 2, 3, 4, 5),
  Observed_Concentration = c(0.072, 2.519, 0.468, 10.951, 11.904),
  Truth_Concentration = c(0, 0.5, 1, 2, 4)
)

# Perform the linear regression for Cysta
regression_model_Cysta <- lm(Observed_Concentration ~ Truth_Concentration, data = data_Cysta)
regression_summary_Cysta <- summary(regression_model_Cysta)

# Plot for Cysta
plot_Cysta <- ggplot(data_Cysta, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
  labs(title = "Actual vs Observed Concentrations of Cysta",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 4) +
  ylim(0, 15) +
  annotate("text", x = 3.5, y = 14, label = paste("R² =", round(regression_summary_Cysta$r.squared, 2), 
                                                  "\nIntercept =", round(regression_model_Cysta$coefficients[1], 2), 
                                                  "\np-value =", round(regression_summary_Cysta$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  theme(panel.border = element_rect(color = "gray", fill = NA, size = 1),
        plot.background = element_rect(color = "gray", fill = NA, size = 1),
        axis.line = element_line(color = "gray"))

# Combine the plots
combined_plot <- grid.arrange(
  plot_2HG, plot_Cysta, 
  ncol = 2, 
  top = textGrob(
    "The Relationship Between Actual and Observed Concentration of 2HG and Cysta (STEAM) at TE = 6ms", 
    gp = gpar(fontsize = 20, fontface = "bold")
  )
)

# Save the combined plot
ggsave("combined_plot_STEAM_6ms.png", plot = combined_plot, width = 12, height = 6)








# concentration analysis at TE 35

# Load necessary library
library(ggplot2)

# Create the data frame with the provided data
data <- data.frame(
  Observed_Concentration = c(5.961, 4.054, 2.3, 10.972, 2.785, 13.05, 6.562, 1.538, 2.335, 9.432, 9.092),
  Actual_Concentration = c(8, 4, 1.6, 12.5, 10, 12.5, 5, 1, 3, 10, 7.5)
)

# Perform the linear regression
regression_model <- lm(Observed_Concentration ~ Actual_Concentration, data = data)
regression_summary <- summary(regression_model)

# Print the regression summary statistics
print(regression_summary)

# Extract the coefficients
intercept <- regression_model$coefficients[1]
slope <- regression_model$coefficients[2]

# Print the regression equation
cat("Regression equation:\n")
cat("Observed Concentration = ", round(intercept, 2), " + ", round(slope, 2), " * Actual Concentration\n")

# Create the plot with the regression line
plot <- ggplot(data, aes(x = Actual_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 2) +
  labs(title = "The Relationship between the  Actual vs Observed Concentrations (2HG 8mM and Cysta 4mM) at PRESS TE = 35ms  ",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  annotate("text", x = max(data$Actual_Concentration), y = max(data$Observed_Concentration), 
           label = paste("R² =", round(regression_summary$r.squared, 2), 
                         "\nIntercept =", round(intercept, 2), 
                         "\nSlope =", round(slope, 2), 
                         "\np-value =", round(regression_summary$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue")

# Display the plot
print(plot)

# Save the plot
ggsave("linear_regression_plot.png", plot = plot, width = 8, height = 6)





# 78

# Creating the data frame with the given data
data_combined <- data.frame(
  Metabolite = c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"),
  Actual_Ratio = c(0.8, 0.4, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75),
  Observed_Ratio = c(0.259, 0.082, 0.104, 1.136, 0, 0.313, 0.235, 0.146, 0.257, 1, 0.512)
)

# Performing linear regression analysis
regression_model <- lm(Observed_Ratio ~ Actual_Ratio, data = data_combined)

# Summary of the regression model
summary(regression_model)

# Extracting regression coefficients
intercept <- coef(regression_model)[1]
slope <- coef(regression_model)[2]

# Display the regression equation
cat("Regression Equation: Observed_Ratio =", round(intercept, 3), "+", round(slope, 3), "* Actual_Ratio\n")

# Plotting the regression line
plot(data_combined$Actual_Ratio, data_combined$Observed_Ratio, 
     main = "Regression Analysis: Actual vs Observed Ratios",
     xlab = "Actual Ratio", ylab = "Observed Ratio",
     pch = 16, col = "blue")
abline(regression_model, col = "red", lwd = 2)






# Correction for RG



# Ensure necessary libraries are loaded
library(dplyr)

# Verification of ATTH2O and WCONC values
ATTH2O <- 1
WCONC <- 55556

if (ATTH2O != 1 | WCONC != 55556) {
  stop("ATTH2O or WCONC values are not appropriate for this phantom study.")
}

# Original concentrations
data <- data.frame(
  Metabolite = c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"),
  Concentration = c(5.961, 4.054, 2.3, 10.972, 2.785, 13.05, 6.562, 1.538, 2.335, 9.432, 9.092)
)

# Correction factor based on Gain Unsuppressed / Gain Suppressed
correction_factor <- 1 / 0.6337  # Inverse to increase concentrations

# Apply the correction factor to the concentrations
data <- data %>%
  mutate(Adjusted_Concentration = Concentration * correction_factor)

# Adjusted Cr concentration
adjusted_cr <- data %>% filter(Metabolite == "Cr") %>% pull(Adjusted_Concentration)

# Calculate new ratios with respect to Cr
data <- data %>%
  mutate(New_Ratio_to_Cr = Adjusted_Concentration / adjusted_cr)

# Print the resulting data frame
print(data)




library(dplyr)

# Constants for ATTH2O and WCONC
ATTH2O <- 1
WCONC <- 55556

# Example Signal Intensity data (this needs to be provided or calculated)
# Placeholder for demonstration purposes
signal_intensity <- c(5.961, 4.054, 2.3, 10.972, 2.785, 13.05, 6.562, 1.538, 2.335, 9.432, 9.092)

# Correction factor based on Gain Unsuppressed / Gain Suppressed
correction_factor <- 1 / 0.6337  # Inverse to increase concentrations

# Calculate concentrations using the provided formula
data <- data.frame(
  Metabolite = c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"),
  Signal_Intensity = signal_intensity
)

data <- data %>%
  mutate(Normalized_Concentration = (Signal_Intensity / (ATTH2O * WCONC)) * correction_factor)

# Print the resulting data frame
print(data)

# Linear regression between Signal Intensity and Normalized Concentration
lm_model <- lm(Normalized_Concentration ~ Signal_Intensity, data = data)
summary(lm_model)

# Generate a regression plot
regression_plot <- ggplot(data, aes(x = Signal_Intensity, y = Normalized_Concentration)) +
  geom_point(size = 5) + 
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
  labs(title = "Regression of Signal Intensity vs Normalized Concentration",
       x = "Signal Intensity",
       y = "Normalized Concentration") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14))

# Print the plot
print(regression_plot)

# Print summary of the linear model
print(summary(lm_model))






















# Create the data frame for 2HG with observed and truth concentrations for PRESS
data_PRESS <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(0, 0.234, 0.488, 0, 1.987),
  Truth_Concentration = c(0, 1, 2, 4, 8)
)

# Perform the linear regression for PRESS
regression_model_PRESS <- lm(Observed_Concentration ~ Truth_Concentration, data = data_PRESS)

# Summary of the regression model to get coefficients and p-values
regression_summary_PRESS <- summary(regression_model_PRESS)

# Extract the p-value for the slope coefficient (Truth_Concentration)
p_value_slope <- regression_summary_PRESS$coefficients[2, 4]
cat("p-value for the slope coefficient:", p_value_slope, "\n")

# Extract t-value for the slope coefficient (Truth_Concentration)
t_value_slope <- regression_summary_PRESS$coefficients[2, 3]
cat("t-value for the slope coefficient:", t_value_slope, "\n")
# Plot the data for PRESS
library(ggplot2)
plot_PRESS <- ggplot(data_PRESS, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "The relationship between the actual and observed concentrations of 2HG (8 mM)  at PRESS TE = 78 ms",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 8) +
  ylim(0, 2) +
  annotate("text", x = 7, y = 1.8, label = paste("R² =", round(regression_summary_PRESS$r.squared, 2), 
                                                 "\nIntercept =", round(regression_model_PRESS$coefficients[1], 2), 
                                                 "\np-value =", round(p_value_slope, 3)),
           size = 5, hjust = 1, vjust = 1, color = "blue")
print(plot_PRESS)
