# best fit analysis and the detection of metabolites using the highest concentration and the PRESS TE = 35, 78, and 97ms 
##### paired t.test for comparison between the true and observed ratio to Cr. 
# checking for false-positive 

library(dplyr)
library(ggplot2)
library(gridExtra)

# Sample data with updated TE observed ratios
data_combined <- data.frame(
  Metabolite = rep(c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"), times = 3),
  Actual_Ratio = rep(c(0.8, 0.4, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75), times = 3),
  Observed_Ratio = c(
    0.603, 0.43, 0.244, 1.163, 0.295, 1.384, 0.696, 0.163, 0.248, 1, 0.964, # TE 35ms
    0.259, 0.082, 0.104, 1.136, 0, 0.313, 0.235, 0.146, 0.257, 1, 0.512, # TE 78ms
    0.205, 0.138, 0.092, 0.911, 0, 0.275, 0.218, 0.102, 0.284, 1, 0.547  # TE 97ms
  ),
  TE = rep(c("35ms", "78ms", "97ms"), each = 11)
)

# Define a function to perform linear regression, summary extraction, and plotting
analyze_te <- function(data, te_label) {
  data_te <- filter(data, TE == te_label)
  lm_model <- lm(Observed_Ratio ~ Actual_Ratio, data = data_te)
  summary_model <- summary(lm_model)
  coef_model <- coef(lm_model)
  r_squared <- summary_model$r.squared
  p_value <- coef(summary_model)[2, 4]
  
  cat(paste("TE =", te_label, ": Observed Ratio =", coef_model[1], "+", coef_model[2], "* Actual Ratio\n"))
  cat(paste("TE =", te_label, ": R² =", round(r_squared, 3), ", p-value =", round(p_value, 3), "\n"))
  
  ggplot(data_te, aes(x = Actual_Ratio, y = Observed_Ratio)) +
    geom_point(size = 7) +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 3) +
    labs(title = paste("True vs Observed Ratio for TE =", te_label), 
         x = "True Ratio", y = "Observed Ratio") +
    annotate("text", x = 1, y = max(data_te$Observed_Ratio) - 0.1, 
             label = paste("R² =", round(r_squared, 3), "\np-value =", round(p_value, 3), 
                           "\ny =", round(coef_model[1], 3), "+", round(coef_model[2], 3), "x"), 
             size = 8, hjust = 0.05, color = "blue") +
    theme_minimal() +
    theme(plot.title = element_text(size = 26, face = "bold"),
          plot.background = element_rect(color = "grey", size = 1),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          panel.grid.major = element_line(size = 1.5, linetype = 'solid', colour = "grey80"),
          panel.grid.minor = element_line(size = 1, linetype = 'solid', colour = "grey90"),
          panel.border = element_rect(colour = "darkgray", fill = NA, size = 2))
}

# Generate and print plots for each TE
plots <- lapply(unique(data_combined$TE), function(te) analyze_te(data_combined, te))
grid.arrange(grobs = plots, ncol = 1)

# Paired T-tests
perform_t_test <- function(te_label) {
  data_te <- filter(data_combined, TE == te_label)
  t_test_result <- t.test(data_te$Actual_Ratio, data_te$Observed_Ratio, paired = TRUE)
  cat(paste("Paired T-test for TE =", te_label, ":\n"))
  print(t_test_result)
}

# Execute paired t-tests for each TE
lapply(unique(data_combined$TE), perform_t_test)


# checking for false-positive 

# dataset TE 35ms 

library(dplyr)

# Create the dataset TE 35ms 
data_combined <- data.frame(
  Metabolite = rep(c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"), times = 3),
  Actual_Ratio = rep(c(0.8, 0.4, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75), times = 3),
  Observed_Ratio = c(
    0.603, 0.43, 0.244, 1.163, 0.295, 1.384, 0.696, 0.163, 0.248, 1, 0.964, # TE 35ms
    0.712, 0.352, 0.177, 1.278, 0.511, 1.295, 0.498, 0.108, 0.272, 0.111, 0.745, # TE 78ms
    0.805, 0.402, 0.162, 1.245, 0.512, 1.238, 0.502, 0.104, 0.305, 0.101, 0.750 # TE 97ms
  ),
  TE = rep(c("35ms", "78ms", "97ms"), each = 11)
)

# Define a threshold for false positives
threshold <- 0.05

# Calculate the difference between Actual_Ratio and Observed_Ratio
data_combined <- data_combined %>%
  mutate(Difference = abs(Actual_Ratio - Observed_Ratio))

# Identify false positives
data_combined <- data_combined %>%
  mutate(False_Positive = ifelse(Difference > threshold, 1, 0))

# Calculate the false positive percentage
false_positive_percentage <- data_combined %>%
  summarise(False_Positive_Percentage = mean(False_Positive) * 100)

# Print the false positive percentage
cat("False Positive Percentage:\n")
print(false_positive_percentage)




# TE=78ms false +


# Create the data frame
data_combined <- data.frame(
  Metabolite = rep(c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"), times = 3),
  Actual_Ratio = rep(c(0.8, 0.4, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75), times = 3),
  Observed_Ratio = c(0.712, 0.352, 0.177, 1.278, 0.511, 1.295, 0.498, 0.108, 0.272, 0.111, 0.745)  # TE 78ms
)

# Define a threshold for false positives
threshold <- 0.05

# Calculate the difference between Actual_Ratio and Observed_Ratio
data_combined <- data_combined %>%
  mutate(Difference = abs(Actual_Ratio - Observed_Ratio))

# Identify false positives
data_combined <- data_combined %>%
  mutate(False_Positive = ifelse(Difference > threshold, 1, 0))

# Calculate the false positive percentage
false_positive_percentage <- data_combined %>%
  summarise(False_Positive_Percentage = mean(False_Positive) * 100)

# Print the false positive percentage
cat("False Positive Percentage:\n")
print(false_positive_percentage)




# TE = 97ms false-postitive 

# Create the data frame
data_combined <- data.frame(
  Metabolite = rep(c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"), times = 1),
  Actual_Ratio = c(0.8, 0.4, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75),
  Observed_Ratio = c(0.205, 0.138, 0.092, 0.911, 0, 0.275, 0.218, 0.102, 0.284, 1, 0.547)  # TE 97ms data
)

# Paired T-test for TE = 97ms
t_test_97 <- t.test(data_combined$Actual_Ratio, data_combined$Observed_Ratio, paired = TRUE)
cat("Paired T-test for TE = 97ms:\n")
print(t_test_97)

# Define a threshold for false positives
threshold <- 0.05

# Calculate the difference between Actual_Ratio and Observed_Ratio
data_combined <- data_combined %>%
  mutate(Difference = abs(Actual_Ratio - Observed_Ratio))

# Identify false positives
data_combined <- data_combined %>%
  mutate(False_Positive = ifelse(Difference > threshold, 1, 0))

# Calculate the false positive percentage
false_positive_percentage <- data_combined %>%
  summarise(False_Positive_Percentage = mean(False_Positive) * 100)

# Print the false positive percentage
cat("False Positive Percentage:\n")
print(false_positive_percentage)


                
# STEAM analysis ofspectral fitting using true ratio and observed ratio to Creatine 


# STEAM analysis of fitting spectra using true ratio and observed ratio TE = 6ms, 4mM
library(dplyr)
library(ggplot2)

# Sample data
data_steam_6ms <- data.frame(
  Metabolite = c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"),
  Actual_Ratio = c(0.4, 0.2, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75),
  Observed_Ratio = c(0.328, 0.334, 0.182, 1.136, 0.248, 1.043, 0.43, 0.241, 0.259, 1, 0.961)
)

# Fit the linear regression model
lm_model_steam_6ms <- lm(Observed_Ratio ~ Actual_Ratio, data = data_steam_6ms)
summary(lm_model_steam_6ms)

# Extract model details for annotation
intercept <- round(coef(lm_model_steam_6ms)[1], 4)
slope <- round(coef(lm_model_steam_6ms)[2], 4)
r_squared <- round(summary(lm_model_steam_6ms)$r.squared * 100, 2)
p_value <- round(summary(lm_model_steam_6ms)$coefficients[2, 4], 5)

# Generate the regression plot
plot_steam_6ms <- ggplot(data_steam_6ms, aes(x = Actual_Ratio, y = Observed_Ratio)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 2) +
  labs(title = "Obsereved vs True Ratio STEAM TE = 6ms (Phantom 2)", x = "True Ratios", y = "Observed Ratios") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        plot.background = element_rect(color = "grey", size = 1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        panel.grid.major = element_line(size = 1.5, linetype = 'solid', colour = "grey80"),
        panel.grid.minor = element_line(size = 1, linetype = 'solid', colour = "grey90"),
        panel.border = element_rect(colour = "darkgray", fill = NA, size = 2)) +
  annotate("text", x = 0.2, y = 1.08, label = paste("y =", intercept, "+", slope, "* x"), size = 4, color = "blue", hjust = -1.3) +
  annotate("text", x = 0.2, y = 1.03, label = paste("R² =", r_squared, "%"), size = 4, color = "blue", hjust = -2.3) +
  annotate("text", x = 0.2, y = 0.99, label = paste("p =", p_value), size = 4, color = "blue", hjust = -2.5)

print(plot_steam_6ms)

# Sample data
true_ratio <- c(0.4, 0.2, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75)
observed_ratio <- c(0.328, 0.334, 0.182, 1.136, 0.248, 1.043, 0.43, 0.241, 0.259, 1, 0.961)

# Perform paired t-test
paired_t_test_result <- t.test(true_ratio, observed_ratio, paired = TRUE)
print(paired_t_test_result)


#false-positive 6ms 4mM


# Define a threshold for false positives
threshold <- 0.05

# Calculate the difference between true_ratio and observed_ratio
difference <- abs(true_ratio - observed_ratio)

# Identify false positives
false_positive <- ifelse(difference > threshold, 1, 0)

# Calculate the false positive percentage
false_positive_percentage <- mean(false_positive) * 100

# Print the false positive percentage
cat("False Positive Percentage for STEAM at TE = 6ms for 4mN:\n")
print(false_positive_percentage)



#####################

# STEAM analysis of fitting spectra using true ratio and observed ratio TE = 6ms and con 8mM
library(dplyr)
library(ggplot2)

# Sample data
data_steam_6ms <- data.frame(
  Metabolite = c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"),
  Actual_Ratio = c(0.8, 0.4, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75),
  Observed_Ratio = c(0.777, 0.557, 0.166, 1.175, 0.327, 1.078, 0.451, 0.306, 0.25, 1, 0.897)
)

# Fit the linear regression model
lm_model_steam_6ms <- lm(Observed_Ratio ~ Actual_Ratio, data = data_steam_6ms)
summary(lm_model_steam_6ms)

# Extract model details for annotation
intercept <- round(coef(lm_model_steam_6ms)[1], 4)
slope <- round(coef(lm_model_steam_6ms)[2], 4)
r_squared <- round(summary(lm_model_steam_6ms)$r.squared * 100, 2)
p_value <- round(summary(lm_model_steam_6ms)$coefficients[2, 4], 5)

# Generate the regression plot
plot_steam_6ms <- ggplot(data_steam_6ms, aes(x = Actual_Ratio, y = Observed_Ratio)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 2) +
  labs(title = "Obsereved vs True Ratio for STEAM TE = 6ms", x = "True Ratios", y = "Observed Ratios") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.background = element_rect(color = "grey", size = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    panel.grid.major = element_line(size = 1.5, linetype = 'solid', colour = "grey80"),
    panel.grid.minor = element_line(size = 1, linetype = 'solid', colour = "grey90"),
    panel.border = element_rect(colour = "darkgray", fill = NA, size = 2)
  ) +
  annotate("text", x = 0.15, y = 1.1, label = paste("y =", intercept, "+", slope, "* x"), size = 4, color = "blue", hjust = -1.3) +
  annotate("text", x = 0.15, y = 1.05, label = paste("R² =", r_squared, "%"), size = 4, color = "blue", hjust = -2.8) +
  annotate("text", x = 0.15, y = 1.00, label = paste("p =", p_value), size = 4, color = "blue", hjust = -3)

print(plot_steam_6ms)


#false-positive 6ms 8mM

# Load necessary library
library(dplyr)

# Define the data
data_steam_6ms <- data.frame(
  Metabolite = c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"),
  Actual_Ratio = c(0.8, 0.4, 0.16, 1.25, 0.5, 1.25, 0.5, 0.1, 0.3, 0.10, 0.75),
  Observed_Ratio = c(0.777, 0.557, 0.166, 1.175, 0.327, 1.078, 0.451, 0.306, 0.25, 1, 0.897)
)

# Define a threshold for false positives
threshold <- 0.05

# Calculate the difference between Actual_Ratio and Observed_Ratio
data_steam_6ms <- data_steam_6ms %>%
  mutate(Difference = abs(Actual_Ratio - Observed_Ratio))

# Identify false positives
data_steam_6ms <- data_steam_6ms %>%
  mutate(False_Positive = ifelse(Difference > threshold, 1, 0))

# Calculate the false positive percentage
false_positive_percentage <- data_steam_6ms %>%
  summarise(False_Positive_Percentage = mean(False_Positive) * 100)

# Print the false positive percentage
cat("False Positive Percentage for STEAM at TE = 6ms:\n")
print(false_positive_percentage)

