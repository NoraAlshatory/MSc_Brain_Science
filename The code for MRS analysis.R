
# best fit analysis and the detection of metabolites using the highest concentration and the PRESS TE = 35, 78, and 97ms 



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

# Fit linear regression models for each TE using direct indexing
lm_model_35 <- lm(Observed_Ratio[1:11] ~ Actual_Ratio[1:11], data = data_combined)
lm_model_78 <- lm(Observed_Ratio[12:22] ~ Actual_Ratio[12:22], data = data_combined)
lm_model_97 <- lm(Observed_Ratio[23:33] ~ Actual_Ratio[23:33], data = data_combined)

# Summarize the models
summary_35 <- summary(lm_model_35)
summary_78 <- summary(lm_model_78)
summary_97 <- summary(lm_model_97)

# Extract coefficients for the equations
coef_35 <- coef(lm_model_35)
coef_78 <- coef(lm_model_78)
coef_97 <- coef(lm_model_97)

# Extract R-squared and p-values
r2_35 <- summary_35$r.squared
p_35 <- coef(summary_35)[2, 4]

r2_78 <- summary_78$r.squared
p_78 <- coef(summary_78)[2, 4]

r2_97 <- summary_97$r.squared
p_97 <- coef(summary_97)[2, 4]

# Print the linear equations
cat("TE = 35ms: Observed Ratio = ", coef_35[1], " + ", coef_35[2], " * Actual Ratio\n")
cat("TE = 78ms: Observed Ratio = ", coef_78[1], " + ", coef_78[2], " * Actual Ratio\n")
cat("TE = 97ms: Observed Ratio = ", coef_97[1], " + ", coef_97[2], " * Actual Ratio\n")

# Print R-squared and p-values
cat("TE = 35ms: R² = ", r2_35, ", p-value = ", p_35, "\n")
cat("TE = 78ms: R² = ", r2_78, ", p-value = ", p_78, "\n")
cat("TE = 97ms: R² = ", r2_97, ", p-value = ", p_97, "\n")

# Generate regression plots with visible grid lines and inside border
plot_35 <- ggplot(data_combined[1:11,], aes(x = as.numeric(as.character(Actual_Ratio)), y = Observed_Ratio)) +
  geom_point(size = 7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 3) +
  labs(title = "True vs Observed Ratio for TE = 35ms", x = "True Ratio", y = "Observed Ratio") +
  annotate("text", x = 1, y = max(data_combined$Observed_Ratio[1:11]) - 0.1, 
           label = paste("R² = ", round(r2_35, 3), "\np-value = ", round(p_35, 3), 
                         "\ny = ", round(coef_35[1], 3), " + ", round(coef_35[2], 3), "x"), 
           size = 8, hjust = 0.05, color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(size = 26, face = "bold"),
        plot.background = element_rect(color = "grey", size = 1),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(size = 1.5, linetype = 'solid', colour = "grey80"),
        panel.grid.minor = element_line(size = 1, linetype = 'solid', colour = "grey90"),
        panel.border = element_rect(colour = "darkgray", fill = NA, size = 2))
print(plot_35)

plot_78 <- ggplot(data_combined[12:22,], aes(x = as.numeric(as.character(Actual_Ratio)), y = Observed_Ratio)) +
  geom_point(size = 7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 3) +
  labs(title = "True vs Observed Ratio for TE = 78ms", x = "True Ratio", y = "Observed Ratio") +
  annotate("text", x = 1, y = max(data_combined$Observed_Ratio[12:22]) - 0.1, 
           label = paste("R² = ", round(r2_78, 3), "\np-value = ", round(p_78, 3), 
                         "\ny = ", round(coef_78[1], 3), " + ", round(coef_78[2], 3), "x"), 
           size = 8, hjust = 0.05, color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(size = 26, face = "bold"),
        plot.background = element_rect(color = "grey", size = 1),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(size = 1.5, linetype = 'solid', colour = "grey80"),
        panel.grid.minor = element_line(size = 1, linetype = 'solid', colour = "grey90"),
        panel.border = element_rect(colour = "darkgray", fill = NA, size = 2))
print(plot_78)


plot_97 <- ggplot(data_combined[23:33,], aes(x = as.numeric(as.character(Actual_Ratio)), y = Observed_Ratio)) +
  geom_point(size = 7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 3) +
  labs(title = "True vs Observed Ratio for TE = 97ms", x = "True Ratio", y = "Observed Ratio") +
  annotate("text", x = 1, y = max(data_combined$Observed_Ratio[23:33]) - 0.1, 
           label = paste("R² = ", round(r2_97, 3), "\np-value = ", round(p_97, 3), 
                         "\ny = ", round(coef_97[1], 3), " + ", round(coef_97[2], 3), "x"), 
           size = 8, hjust = 0.05, color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(size = 26, face = "bold"),
        plot.background = element_rect(color = "grey", size = 1),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(size = 1.5, linetype = 'solid', colour = "grey80"),
        panel.grid.minor = element_line(size = 1, linetype = 'solid', colour = "grey90"),
        panel.border = element_rect(colour = "darkgray", fill = NA, size = 2))
print(plot_97)




##### paired t.test for comparison between the true and observed ratio at different TE for phantom 1

# Create the data frame
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

# Paired T-test for TE = 35ms
data_TE_35 <- subset(data_combined, TE == "35ms")
t_test_35 <- t.test(data_TE_35$Actual_Ratio, data_TE_35$Observed_Ratio, paired = TRUE)
cat("Paired T-test for TE = 35ms:\n")
print(t_test_35)

# Paired T-test for TE = 78ms
data_TE_78 <- subset(data_combined, TE == "78ms")
t_test_78 <- t.test(data_TE_78$Actual_Ratio, data_TE_78$Observed_Ratio, paired = TRUE)
cat("Paired T-test for TE = 78ms:\n")
print(t_test_78)

# Paired T-test for TE = 97ms
data_TE_97 <- subset(data_combined, TE == "97ms")
t_test_97 <- t.test(data_TE_97$Actual_Ratio, data_TE_97$Observed_Ratio, paired = TRUE)
cat("Paired T-test for TE = 97ms:\n")
print(t_test_97)






#cheaking for false postive 

# Load necessary libraries
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


# Load necessary libraries
library(dplyr)

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






# TE = 97ms false +
# Load necessary libraries
library(dplyr)

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


###################

# STEAM analysis of fitting spectra using true ratio and observed ratio 


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


#false + 6ms 4mM


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
  labs(title = "Obsereved vs True Ratio for STEAM TE = 6ms (Phantom 1)", x = "True Ratios", y = "Observed Ratios") +
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




#false + 6ms 8mM

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




#### bland Altman plot of the observed ratios across all TE = 35, 78, and 97ms for PRERSS and STEAM



##### for PRESS

library(dplyr)
library(ggplot2)

# Sample data 
data <- data.frame(
  Metabolite = rep(c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"), 15),
  Condition = rep(c("PRESS_0mM_35", "PRESS_1mM_35", "PRESS_2mM_35", "PRESS_4mM_35", "PRESS_8mM_35",
                    "PRESS_0mM_78", "PRESS_1mM_78", "PRESS_2mM_78", "PRESS_4mM_78", "PRESS_8mM_78",
                    "PRESS_0mM_97", "PRESS_1mM_97", "PRESS_2mM_97", "PRESS_4mM_97", "PRESS_8mM_97"), each = 11),
  True_Ratio = c(0, 0, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.1, 0.05, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.2, 0.1, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.4, 0.2, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.8, 0.4, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0, 0, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.1, 0.05, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.2, 0.1, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.4, 0.2, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.8, 0.4, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0, 0, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.1, 0.05, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.2, 0.1, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.4, 0.2, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.8, 0.4, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75), 
  Observed_Ratio = c(0.525, 0.903, 0.253, 0.903, 0.336, 1.185, 0.417, 0.0015, 0.275, 1, 1.266,
                     1.836, 0.158, 0.254, 1.017, 0.384, 1.474, 0.668, 0.012, 0.255, 1, 1.022,
                     0, 0.155, 0.219, 0.967, 0.324, 1.135, 0.498, 0, 0.253, 1, 1.152,
                     0.963, 0.172, 0.197, 0.797, 0.285, 0.757, 0.328, 0.0, 0.249, 1, 1.097,
                     0.603, 0.43, 0.244, 1.163, 0.295, 1.384, 0.696, 0.163, 0.248, 1, 0.964,
                     0, 0.065, 0.125, 0.894, 0, 0.241, 0.235, 0.212, 0.298, 1, 0.616,
                     0.099, 0.027, 0.093, 1.175, 0, 0.358, 0.246, 0.191, 0.281, 1, 0.643,
                     0.077, 0.082, 0.141, 0.962, 0, 0.277, 0.253, 0.153, 0.282, 1, 0.564,
                     0, 0.122, 0.131, 0.906, 0, 0.253, 0.253, 0.193, 0.291, 1, 0.627,
                     0.259, 0.082, 0.104, 1.136, 0, 0.313, 0.235, 0.146, 0.257, 1, 0.512,
                     0.052, 0.158, 0.165, 0.897, 0, 0.305, 0.246, 0.126, 0.305, 1, 0.594,
                     0.011, 0.087, 0.056, 0.832, 0, 0.278, 0.201, 0.132, 0.277, 1, 0.578,
                     0.037, 0.093, 0.141, 0.834, 0, 0.273, 0.241, 0.091, 0.29, 1, 0.592,
                     0.136, 0.188, 0.153, 0.907, 0, 0.335, 0.272, 0.128, 0.297, 1, 0.589,
                     0.205, 0.138, 0.092, 0.911, 0, 0.275, 0.218, 0.102, 0.284, 1, 0.547)
)

# Calculate the mean and difference between true and observed ratios
data <- data %>%
  mutate(Mean = (Observed_Ratio + True_Ratio) / 2,
         Difference = Observed_Ratio - True_Ratio)

# Calculate the mean difference and the limits of agreement
mean_diff <- mean(data$Difference)
sd_diff <- sd(data$Difference)
upper_limit <- mean_diff + 1.96 * sd_diff
lower_limit <- mean_diff - 1.96 * sd_diff

# Create the Bland-Altman plot
plot <- ggplot(data, aes(x = Mean, y = Difference)) +
  geom_point(color = "black", size = 3.5) +
  geom_hline(yintercept = mean_diff, color = "blue", linetype = "dashed", size = 2) +
  geom_hline(yintercept = upper_limit, color = "red", linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = lower_limit, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Bland-Altman Plot of the Observed Ratio vs True Ratio of PRESS Sequence",
       x = "Mean",
       y = "Difference") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(face = "bold"),
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1)
  ) +
  annotate("text", x = Inf, y = mean_diff, label = paste("Mean:", round(mean_diff, 2)), vjust = -1, hjust = 1.1, size = 6, color = "red") +
  annotate("text", x = Inf, y = upper_limit, label = paste("Upper Limit:", round(upper_limit, 2)), vjust = -1, hjust = 1.1, size = 6, color = "black") +
  annotate("text", x = Inf, y = lower_limit, label = paste("Lower Limit:", round(lower_limit, 2)), vjust = 1.1, hjust = 1.1, size = 6, color = "black")

# Print the plot to display it
print(plot)

# Perform a paired t-test
t_test_result <- t.test(data$Observed_Ratio, data$True_Ratio, paired = TRUE)

# Print the results
print(t_test_result)





# bland Altman plot for STEAM 



# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Sample data
data <- data.frame(
  Metabolite = rep(c("2HG", "Cysta", "NAAG", "NAA", "Gln", "Glu", "Lac", "Gly", "Cho", "Cr", "ml"), 15),
  Condition = rep(c("PRESS_0mM_35", "PRESS_1mM_35", "PRESS_2mM_35", "PRESS_4mM_35", "PRESS_8mM_35",
                    "PRESS_0mM_78", "PRESS_1mM_78", "PRESS_2mM_78", "PRESS_4mM_78", "PRESS_8mM_78",
                    "PRESS_0mM_97", "PRESS_1mM_97", "PRESS_2mM_97", "PRESS_4mM_97", "PRESS_8mM_97"), each = 11),
  True_Ratio = c(0, 0, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.1, 0.05, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.2, 0.1, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.4, 0.2, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.8, 0.4, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0, 0.0, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.1, 0.05, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.2, 0.1, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.4, 0.2, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.8, 0.4, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0, 0, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.1, 0.05, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.2, 0.1, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.4, 0.2, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75, 
                 0.8, 0.4, 0.16, 0.125, 0.5, 0.125, 0.5, 0.1, 0.3, 0.1, 0.75),
  Observed_Ratio = c(0.509, 0.391, 0.217, 1.042, 0.56, 1.385, 0.305, 0.197, 0.278, 1, 1.347,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0.521, 0.35, 0.218, 1.036, 0.552, 1.31, 0.317, 0.164, 0.273, 1, 1.306,
                     0.775, 0.51, 0.214, 1.027, 0.492, 1.338, 0.31, 0.181, 0.27, 1, 1.249,
                     1.01, 0.681, 0.276, 1.075, 0.538, 1.271, 0.321, 0.213, 0.267, 1, 1.215,
                     2.068, 0, 0.196, 0.985, 0.012, 0, 0.154, 0, 0.291, 1, 1.17,
                     2.197, 0, 0.216, 1.055, 0.254, 0.26, 0.268, 0, 0.282, 1, 1.06,
                     2.496, 0.956, 0.099, 0.924, 0.082, 0.024, 0.156, 0, 0.256, 1, 0.978,
                     2.08, 0, 0.197, 0.98, 0.147, 0.065, 0.178, 0, 0.284, 1, 1.092,
                     2.928, 1.164, 0.14, 0.978, 0.08, 0, 0.149, 0, 0.252, 1, 0.914,
                     0.571, 0.582, 0.194, 0.924, 0, 0.265, 0.066, 0, 0.286, 1, 0.993,
                     0, 0, 0.211, 1.044, 0, 0.541, 0.231, 0, 0.27, 1, 0.973,
                     3.061, 0.46, 0.31, 0.974, 0.01, 0.949, 0.068, 0.04, 0.287, 1, 0.793,
                     2.659, 0.655, 0.261, 0.948, 0.05, 0.73, 0.099, 0.059, 0.279, 1, 0.773,
                     3.234, 0.517, 0.323, 1.021, 0, 1.215, 0.103, 0.057, 0.28, 1, 0.755)
)

# Calculate the mean and difference between true and observed ratios
data <- data %>%
  mutate(Mean = (Observed_Ratio + True_Ratio) / 2,
         Difference = Observed_Ratio - True_Ratio)

# Calculate the mean difference and the limits of agreement
mean_diff <- mean(data$Difference)
sd_diff <- sd(data$Difference)
upper_limit <- mean_diff + 1.96 * sd_diff
lower_limit <- mean_diff - 1.96 * sd_diff

# Create the Bland-Altman plot
plot <- ggplot(data, aes(x = Mean, y = Difference)) +
  geom_point(color = "black", size = 3.5) +
  geom_hline(yintercept = mean_diff, color = "blue", linetype = "dashed", size = 2) +
  geom_hline(yintercept = upper_limit, color = "red", linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = lower_limit, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Bland-Altman Plot of the Observed Ratio vs True Ratio of STEAM Sequence",
       x = "Mean",
       y = "Difference") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(face = "bold"),
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1)
  ) +
  annotate("text", x = Inf, y = mean_diff, label = paste("Mean:", round(mean_diff, 2)), vjust = -1, hjust = 1.1, size = 6, color = "red") +
  annotate("text", x = Inf, y = upper_limit, label = paste("Upper Limit:", round(upper_limit, 2)), vjust = -1, hjust = 1.1, size = 6, color = "black") +
  annotate("text", x = Inf, y = lower_limit, label = paste("Lower Limit:", round(lower_limit, 2)), vjust = 1.1, hjust = 1.1, size = 6, color = "black")

# Print the plot to display it
print(plot)



# Check for normality of the differences
shapiro_test <- shapiro.test(data$Difference)
print(shapiro_test)


t_test <- t.test(data$Observed_Ratio, data$True_Ratio, paired = TRUE)
print(t_test)





#### regression analysis of the observed and actual concentration across all TEs at different concentrations.

#2HG and Cysta at 6ms 


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
combined_plot <- grid.arrange(plot_2HG, plot_Cysta, ncol = 2, top = "The Relationship Between Actual and Observed Concentration of 2HG and Cysta (STEAM) at TE = 6ms")

# Save the combined plot
ggsave("combined_plot_STEAM_6ms.png", plot = combined_plot, width = 12, height = 6)







## regression analysis of 2HG at  TE = 78ms concentration 


# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)  # This library provides the textGrob function

# Create the data frame for 2HG with observed and truth concentrations for PRESS
data_PRESS <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(0, 0.234, 0.488, 0, 1.987),
  Truth_Concentration = c(0, 1, 2, 4, 8)
)

# Create the data frame for 2HG with observed and truth concentrations for STEAM
data_STEAM <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(29.951, 13.583, 48.174, 51.997, 56.141),
  Truth_Concentration = c(0, 1, 2, 4, 8)
)

# Perform the linear regression for PRESS
regression_model_PRESS <- lm(Observed_Concentration ~ Truth_Concentration, data = data_PRESS)
regression_summary_PRESS <- summary(regression_model_PRESS)

# Perform the linear regression for STEAM
regression_model_STEAM <- lm(Observed_Concentration ~ Truth_Concentration, data = data_STEAM)
regression_summary_STEAM <- summary(regression_model_STEAM)

# Define a custom theme for gray borders
custom_theme <- theme(
  panel.border = element_rect(color = "gray", fill = NA, size = 1),
  plot.background = element_rect(color = "gray", fill = NA, size = 1),
  axis.line = element_line(color = "gray")
)

# Plot the data for PRESS
plot_PRESS <- ggplot(data_PRESS, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "PRESS",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 8) +
  ylim(0, 2) +
  annotate("text", x = 7, y = 1.8, label = paste("R² =", round(regression_summary_PRESS$r.squared, 2), 
                                                 "\nIntercept =", round(regression_model_PRESS$coefficients[1], 2), 
                                                 "\np-value =", round(regression_summary_PRESS$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Plot the data for STEAM
plot_STEAM <- ggplot(data_STEAM, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "STEAM",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 15) +
  xlim(0, 8) +
  ylim(0, 60) +
  annotate("text", x = 7, y = 45, label = paste("R² =", round(regression_summary_STEAM$r.squared, 2), 
                                                "\nIntercept =", round(regression_model_STEAM$coefficients[1], 2), 
                                                "\np-value =", round(regression_summary_STEAM$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Combine the plots with an overarching title
combined_plot <- grid.arrange(
  plot_PRESS, plot_STEAM, 
  ncol = 2, 
  top = textGrob(" (A)  Actual vs Observed Concentrations of 2HG for PRESS and STEAM at TE = 78ms", gp = gpar(fontsize = 20, fontface = "bold"))
)

# Save the combined plot
ggsave("combined_plot.png", plot = combined_plot, width = 12, height = 6)





# cysta at TE = 78ms 
# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)  # This library provides the textGrob function

# Create the data frame for Cysta with observed and truth concentrations for PRESS
data_PRESS <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(1.102, 0.0642, 0.559, 3.87, 0.659),
  Truth_Concentration = c(0, 0.5, 1, 2, 4)
)

# Create the data frame for Cysta with observed and truth concentrations for STEAM
data_STEAM <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(0, 0, 18.452, 0, 22.327),
  Truth_Concentration = c(0, 0.5, 1, 2, 4)
)

# Perform the linear regression for PRESS
regression_model_PRESS <- lm(Observed_Concentration ~ Truth_Concentration, data = data_PRESS)
regression_summary_PRESS <- summary(regression_model_PRESS)

# Perform the linear regression for STEAM
regression_model_STEAM <- lm(Observed_Concentration ~ Truth_Concentration, data = data_STEAM)
regression_summary_STEAM <- summary(regression_model_STEAM)

# Define a custom theme for gray borders
custom_theme <- theme(
  panel.border = element_rect(color = "gray", fill = NA, size = 1),
  plot.background = element_rect(color = "gray", fill = NA, size = 1),
  axis.line = element_line(color = "gray")
)

# Plot the data for PRESS
plot_PRESS <- ggplot(data_PRESS, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "PRESS",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 4) +
  ylim(0, 4) +
  annotate("text", x = 3.5, y = 3.5, label = paste("R² =", round(regression_summary_PRESS$r.squared, 2), 
                                                   "\nIntercept =", round(regression_model_PRESS$coefficients[1], 2), 
                                                   "\np-value =", round(regression_summary_PRESS$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Plot the data for STEAM
plot_STEAM <- ggplot(data_STEAM, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "STEAM",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 15) +
  xlim(0, 4) +
  ylim(0, 24) +
  annotate("text", x = 3.5, y = 20, label = paste("R² =", round(regression_summary_STEAM$r.squared, 2), 
                                                  "\nIntercept =", round(regression_model_STEAM$coefficients[1], 2), 
                                                  "\np-value =", round(regression_summary_STEAM$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Combine the plots with an overarching title
combined_plot <- grid.arrange(
  plot_PRESS, plot_STEAM, 
  ncol = 2, 
  top = textGrob(" (B)  Actual vs Observed Concentrations of Cysta for PRESS and STEAM at TE = 78ms", gp = gpar(fontsize = 20, fontface = "bold"))
)

# Save the combined plot
ggsave("combined_plot.png", plot = combined_plot, width = 12, height = 6)






# TE = 97ms 

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)  # This library provides the textGrob function

# Create the data frame for 2HG with observed and truth concentrations for PRESS
data_PRESS <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(0.792, 0.138, 0.902, 3.962, 5.211),
  Truth_Concentration = c(0, 1, 2, 4, 8)
)

# Create the data frame for 2HG with observed and truth concentrations for STEAM
data_STEAM <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(7.199, 0, 45.472, 58.345, 45.519),
  Truth_Concentration = c(0, 1, 2, 4, 8)
)

# Perform the linear regression for PRESS
regression_model_PRESS <- lm(Observed_Concentration ~ Truth_Concentration, data = data_PRESS)
regression_summary_PRESS <- summary(regression_model_PRESS)

# Perform the linear regression for STEAM
regression_model_STEAM <- lm(Observed_Concentration ~ Truth_Concentration, data = data_STEAM)
regression_summary_STEAM <- summary(regression_model_STEAM)

# Define a custom theme for gray borders
custom_theme <- theme(
  panel.border = element_rect(color = "gray", fill = NA, size = 1),
  plot.background = element_rect(color = "gray", fill = NA, size = 1),
  axis.line = element_line(color = "gray")
)

# Plot the data for PRESS
plot_PRESS <- ggplot(data_PRESS, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "PRESS",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 8) +
  ylim(0, 8) +
  annotate("text", x = 7, y = 7, label = paste("R² =", round(regression_summary_PRESS$r.squared, 2), 
                                               "\nIntercept =", round(regression_model_PRESS$coefficients[1], 2), 
                                               "\np-value =", round(regression_summary_PRESS$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Plot the data for STEAM
plot_STEAM <- ggplot(data_STEAM, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "STEAM",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 8) +
  ylim(0, 60) +
  annotate("text", x = 7, y = 45, label = paste("R² =", round(regression_summary_STEAM$r.squared, 2), 
                                                "\nIntercept =", round(regression_model_STEAM$coefficients[1], 2), 
                                                "\np-value =", round(regression_summary_STEAM$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Combine the plots with an overarching title
combined_plot <- grid.arrange(
  plot_PRESS, plot_STEAM, 
  ncol = 2, 
  top = textGrob("(C)  Actual vs Observed Concentrations of 2HG for PRESS and STEAM at TE = 97ms", gp = gpar(fontsize = 20, fontface = "bold"))
)

# Save the combined plot
ggsave("combined_plot.png", plot = combined_plot, width = 12, height = 6)



# cysta TE =97ms 


# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)  # This library provides the textGrob function

# Create the data frame for Cysta with observed and truth concentrations for PRESS
data_PRESS <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(2.397, 1.089, 2.247, 5.5, 3.511),
  Truth_Concentration = c(0, 0.5, 1, 2, 4)
)

# Create the data frame for Cysta with observed and truth concentrations for STEAM
data_STEAM <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(7.327, 0, 6.831, 14.369, 7.273),
  Truth_Concentration = c(0, 0.5, 1, 2, 4)
)

# Perform the linear regression for PRESS
regression_model_PRESS <- lm(Observed_Concentration ~ Truth_Concentration, data = data_PRESS)
regression_summary_PRESS <- summary(regression_model_PRESS)

# Perform the linear regression for STEAM
regression_model_STEAM <- lm(Observed_Concentration ~ Truth_Concentration, data = data_STEAM)
regression_summary_STEAM <- summary(regression_model_STEAM)

# Define a custom theme for gray borders
custom_theme <- theme(
  panel.border = element_rect(color = "gray", fill = NA, size = 1),
  plot.background = element_rect(color = "gray", fill = NA, size = 1),
  axis.line = element_line(color = "gray")
)

# Plot the data for PRESS
plot_PRESS <- ggplot(data_PRESS, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
  labs(title = "PRESS",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 4) +
  ylim(0, 8) +
  annotate("text", x = 3.5, y = 7.5, label = paste("R² =", round(regression_summary_PRESS$r.squared, 2), 
                                                   "\nIntercept =", round(regression_model_PRESS$coefficients[1], 2), 
                                                   "\nSlope =", round(regression_model_PRESS$coefficients[2], 2), 
                                                   "\np-value =", round(regression_summary_PRESS$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Plot the data for STEAM
plot_STEAM <- ggplot(data_STEAM, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
  labs(title = "STEAM",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 4) +
  ylim(0, 16) +
  annotate("text", x = 3.5, y = 15, label = paste("R² =", round(regression_summary_STEAM$r.squared, 2), 
                                                  "\nIntercept =", round(regression_model_STEAM$coefficients[1], 2), 
                                                  "\nSlope =", round(regression_model_STEAM$coefficients[2], 2), 
                                                  "\np-value =", round(regression_summary_STEAM$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Combine the plots with an overarching title
combined_plot <- grid.arrange(
  plot_PRESS, plot_STEAM, 
  ncol = 2, 
  top = textGrob("(D) Actual vs Observed Concentrations of Cysta for PRESS and STEAM at TE = 97ms", gp = gpar(fontsize = 20, fontface = "bold"))
)

# Save the combined plot
ggsave("combined_plot.png", plot = combined_plot, width = 12, height = 6)

# Print the regression summary statistics
cat("Regression summary for PRESS:\n")
print(regression_summary_PRESS)

cat("\nRegression summary for STEAM:\n")
print(regression_summary_STEAM)

# Print the regression equations
cat("Regression equation for PRESS:\n")
cat("Observed Concentration = ", round(regression_model_PRESS$coefficients[1], 2), " + ", round(regression_model_PRESS$coefficients[2], 2), " * Actual Concentration\n\n")

cat("Regression equation for STEAM:\n")
cat("Observed Concentration = ", round(regression_model_STEAM$coefficients[1], 2), " + ", round(regression_model_STEAM$coefficients[2], 2), " * Actual Concentration\n")





# 144 and 288ms

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)  # This library provides the textGrob function

# Create the data frame for 2HG with observed and truth concentrations for PRESS at 144ms
data_PRESS_144ms <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(13.62, 0, 17.396, 20.012, 17.858),
  Truth_Concentration = c(0, 1, 2, 4, 8)
)

# Create the data frame for 2HG with observed and truth concentrations for PRESS at 288ms
data_PRESS_288ms <- data.frame(
  Phantom = c(0, 4, 3, 2, 1),
  Observed_Concentration = c(1.563, 0, 0, 3.783, 0),
  Truth_Concentration = c(0, 1, 2, 4, 8)
)

# Perform the linear regression for PRESS at 144ms
regression_model_PRESS_144ms <- lm(Observed_Concentration ~ Truth_Concentration, data = data_PRESS_144ms)
regression_summary_PRESS_144ms <- summary(regression_model_PRESS_144ms)

# Perform the linear regression for PRESS at 288ms
regression_model_PRESS_288ms <- lm(Observed_Concentration ~ Truth_Concentration, data = data_PRESS_288ms)
regression_summary_PRESS_288ms <- summary(regression_model_PRESS_288ms)

# Define a custom theme for gray borders
custom_theme <- theme(
  panel.border = element_rect(color = "gray", fill = NA, size = 1),
  plot.background = element_rect(color = "gray", fill = NA, size = 1),
  axis.line = element_line(color = "gray")
)

# Plot the data for PRESS at 144ms
plot_PRESS_144ms <- ggplot(data_PRESS_144ms, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "PRESS at TE = 144ms",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 8) +
  ylim(0, 25) +
  annotate("text", x = 7, y = 22, label = paste("R² =", round(regression_summary_PRESS_144ms$r.squared, 2), 
                                                "\nIntercept =", round(regression_model_PRESS_144ms$coefficients[1], 2), 
                                                "\nSlope =", round(regression_model_PRESS_144ms$coefficients[2], 2), 
                                                "\np-value =", round(regression_summary_PRESS_144ms$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Plot the data for PRESS at 288ms
plot_PRESS_288ms <- ggplot(data_PRESS_288ms, aes(x = Truth_Concentration, y = Observed_Concentration)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, fullrange = TRUE) +
  labs(title = "PRESS at TE = 288ms",
       x = "Actual Concentration (mM)",
       y = "Observed Concentration (mM)") +
  theme_minimal(base_size = 14) +
  xlim(0, 8) +
  ylim(0, 5) +
  annotate("text", x = 7, y = 4.5, label = paste("R² =", round(regression_summary_PRESS_288ms$r.squared, 2), 
                                                 "\nIntercept =", round(regression_model_PRESS_288ms$coefficients[1], 2), 
                                                 "\nSlope =", round(regression_model_PRESS_288ms$coefficients[2], 2), 
                                                 "\np-value =", round(regression_summary_PRESS_288ms$coefficients[2, 4], 2)),
           size = 5, hjust = 1, vjust = 1, color = "blue") +
  custom_theme

# Combine the plots with an overarching title
combined_plot <- grid.arrange(
  plot_PRESS_144ms, plot_PRESS_288ms, 
  ncol = 2, 
  top = textGrob("Regression Analysis Between Actual and Observed Concentrations of 2HG for PRESS at TE = 144ms and 288ms", gp = gpar(fontsize = 20, fontface = "bold"))
)

# Save the combined plot
ggsave("combined_plot_144ms_288ms.png", plot = combined_plot, width = 12, height = 6)

# Print the regression equations
cat("Regression equation for PRESS at 144ms:\n")
cat("Observed Concentration = ", round(regression_model_PRESS_144ms$coefficients[1], 2), " + ", round(regression_model_PRESS_144ms$coefficients[2], 2), " * Actual Concentration\n\n")

cat("Regression equation for PRESS at 288ms:\n")
cat("Observed Concentration = ", round(regression_model_PRESS_288ms$coefficients[1], 2), " + ", round(regression_model_PRESS_288ms$coefficients[2], 2), " * Actual Concentration\n")





# finding the specificity and sensitivity between PRESS and STEAM at TE =35, 78 and 97ms.
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




################################################################

# spectra

###############################################################




# phantom 0 at all TE for STEAM



#phantom 0 at all TE for PRESS  0

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom0_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 3000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.8) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25),limits = c(5, 1)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Control Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")





# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom0_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 3000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.8) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25), limits = c(5, 1)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Control Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")








# spectra close up from 2-3ppm of PRESS Phantom 0

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom0_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation and thicker lines, focusing on x-axis from 2 to 3
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 1.2) +  # Increase the line size
  scale_x_reverse(limits = c(3, 2, by = 0.25), breaks = seq(2, 3, by = 0.25)) +  # Set specific limits and breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")








# phantom 0 at all TE for STEAM

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom0_STEAM.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "6" ~ 6 * 4000,
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("97", "78", "35", "6")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.9) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25), limits = c(5, 1) ) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Control Spectra at Different TEs (STEAM)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 28, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 19, face = "bold"),
        legend.text = element_text(size = 19, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")





# sTEAM at 0 for 2- 3 spectra 



# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom0_STEAM.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 2.00 and 3.00
data_long <- data_long %>%
  filter(PPM >= 2.00 & PPM <= 3.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "6" ~ 6 * 4000,
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("97", "78", "35", "6")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 1.2) +
  scale_x_reverse(breaks = seq(2, 3, by = 0.25), limits = c(3, 2, by = 0.25)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Control Spectra at Different TEs (STEAM)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 1, face = "bold"),
        axis.title = element_text(size = 27, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 27, face = "bold"),
        axis.text.y = element_text(size = 27, face = "bold"),
        legend.title = element_text(size = 27, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")



#####################################################################################################
#####################################################################################################



#





#phantom 1 at all TE for PRESS  1




# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom1_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 3000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.8) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25), limits = c(5, 1)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 1 Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")








# spectra close up from 2-3ppm of PRESS Phantom 1

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom1_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation and thicker lines, focusing on x-axis from 2 to 3
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 1.2) +  # Increase the line size
  scale_x_reverse(limits = c(3, 2, by = 0.25), breaks = seq(2, 3, by = 0.25)) +  # Set specific limits and breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 1 Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")








# phantom 1 at all TE for STEAM

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom1_STEAM.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "6" ~ 6 * 4000,
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("97", "78", "35", "6")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.5) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25), limits = c(5, 1) ) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 1 Spectra at Different TEs (STEAM)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 28, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 19, face = "bold"),
        legend.text = element_text(size = 19, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")





# sTEAM at 1 for 2- 3 spectra 



# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom1_STEAM.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 2.00 and 3.00
data_long <- data_long %>%
  filter(PPM >= 2.00 & PPM <= 3.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "6" ~ 6 * 4000,
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("97", "78", "35", "6")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.5) +
  scale_x_reverse(breaks = seq(2, 3, by = 0.25), limits = c(3, 2, by = 0.25)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 1 Spectra at Different TEs (STEAM)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 1, face = "bold"),
        axis.title = element_text(size = 27, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 27, face = "bold"),
        axis.text.y = element_text(size = 27, face = "bold"),
        legend.title = element_text(size = 27, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")





######################################################################################################

#phantom 2 at all TE for PRESS  2

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom2_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 6000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.9) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 2 Separated Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")





# spectra close up from 2-3ppm of PRESS Phantom 2 

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom2_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation and thicker lines, focusing on x-axis from 2 to 3
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 1) +  # Increase the line size
  scale_x_reverse(limits = c(3, 2), breaks = seq(2, 3, by = 0.25)) +  # Set specific limits and breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 2 Separated Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")







# phantom 2 at all TE for STEAM


library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom2_STEAM.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.8) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 2 Separated Spectra at Different TEs (STEAM)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")












##############################################################################################






#phantom 3 at all TE for PRESS  3

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom3_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.8) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 3 Separated Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")





# spectra close up from 2-3ppm of PRESS Phantom 3

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom3_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation and thicker lines, focusing on x-axis from 2 to 3
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.6) +  # Increase the line size
  scale_x_reverse(limits = c(3, 2), breaks = seq(2, 3, by = 0.1)) +  # Set specific limits and breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 3 Separated Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")






# phantom 3 at all TE for STEAM


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom3_STEAM.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 2.00 and 3.00
data_long <- data_long %>%
  filter(PPM >= 2.00 & PPM <= 3.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "6" ~ 6 * 4000,
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("97", "78", "35", "6")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.5) +
  scale_x_reverse(breaks = seq(2, 3, by = 0.1), limits = c(3, 2)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 3 Separated Spectra at Different TEs (STEAM)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30, face = "bold"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")





########################################################################################






#phantom 4 at all TE for PRESS  4

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom4_PRESS.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "35" ~ 35 * 5500,
    TE == "78" ~ 78 * 4000,
    TE == "97" ~ 97 * 4000,  # Larger offset for TE 97
    TE == "144" ~ 144 * 4000, # Larger offset for TE 144
    TE == "288" ~ 288 * 2500,
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("288", "144", "97", "78", "35")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.8) +
  scale_x_reverse(breaks = seq(1, 7, by = 0.25)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 4 Separated Spectra at Different TEs (PRESS)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")





# spectra close up from 2-3ppm of PRESS Phantom 4

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Phantom4_STEAM.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "TE"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 0.50 and 5.00
data_long <- data_long %>%
  filter(PPM >= 0.50 & PPM <= 5.00)

# Adjust the vertical spacing by adding an offset based on TE
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    TE == "6" ~ 6 * 4000,
    TE == "35" ~ 35 * 4000,
    TE == "78" ~ 78 * 4000,  # Larger offset for TE 97
    TE == "97" ~ 97 * 4000, # Larger offset for TE 144
    TRUE ~ as.numeric(TE) * 4000  # Default case
  ))

# Define the order of TEs for the legend
te_order <- c("97", "78","35", "6")

# Plot the spectra with separation and thicker lines, focusing on x-axis from 2 to 3
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.8) +  # Increase the line size
  scale_x_reverse(limits = c(3, 2), breaks = seq(2, 3, by = 0.1)) +  # Set specific limits and breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Phantom 4 Separated Spectra at Different TEs (STEAM)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")







###########concentration spectra

########################################################################################


#  78 spectra
#########################################################



# spectra concentration at 78 for all concentraios but not at 1mM as its bad

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Con-spectra_data_78.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "Concentration"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 1.00 and 5.00 and remove concentration 1
data_long <- data_long %>%
  filter(PPM >= 1.00 & PPM <= 5.00 & Concentration != "1")

# Adjust the vertical spacing by adding an offset based on Concentration
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    Concentration == "8" ~ 8 * 80000,
    Concentration == "4" ~ 4 * 80000,
    Concentration == "2" ~ 2 * 90000,
    Concentration == "0" ~ 0 * 400,
    TRUE ~ as.numeric(Concentration) * 100000  # Default case
  ))

# Define the order of Concentrations for the legend
concentration_order <- c("8", "4", "2", "0")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(Concentration, levels = concentration_order))) +
  geom_line(size = 0.5) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25), limits = c(5, 1)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Adjust y-axis limits
  labs(title = "The Effect of Concentration on the Spectra TE78ms ",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "Concentration (mM)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")



# from the 2- 3 spectra
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Con-spectra_data_78.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "Concentration"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 2.00 and 3.00 and remove concentration 1
data_long <- data_long %>%
  filter(PPM >= 2.00 & PPM <= 3.00 & Concentration != "1")

# Adjust the vertical spacing by adding an offset based on Concentration
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    Concentration == "8" ~ 8 * 80000,
    Concentration == "4" ~ 4 * 80000,
    Concentration == "2" ~ 2 * 90000,
    Concentration == "0" ~ 0 * 400,
    TRUE ~ as.numeric(Concentration) * 100000  # Default case
  ))

# Define the order of Concentrations for the legend
concentration_order <- c("8", "4", "2", "0")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(Concentration, levels = concentration_order))) +
  geom_line(size = 0.8) +
  scale_x_reverse(breaks = seq(2, 3, by = 0.1), limits = c(3, 2)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Adjust y-axis limits
  labs(title = "The Effect of Concentration on the Spectra TE 78ms",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "Concentration (mM)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot_2_3_PPM.png")







# spectra concentration at 97 for all concentraios but not at 1mM as its bad

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Con_spectra_data_97.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "Concentration"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 1.00 and 5.00 and remove concentration 1
data_long <- data_long %>%
  filter(PPM >= 1.00 & PPM <= 5.00 & Concentration != "1")

# Adjust the vertical spacing by adding an offset based on Concentration
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    Concentration == "8" ~ 8 * 80000,
    Concentration == "4" ~ 4 * 80000,
    Concentration == "2" ~ 2 * 90000,
    Concentration == "0" ~ 0 * 400,
    TRUE ~ as.numeric(Concentration) * 100000  # Default case
  ))

# Define the order of Concentrations for the legend
concentration_order <- c("8", "4", "2", "0")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(Concentration, levels = concentration_order))) +
  geom_line(size = 0.5) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25), limits = c(5, 1)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Adjust y-axis limits
  labs(title = "The Effect of Concentration on the Spectra TE97ms ",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "Concentration (mM)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")



# from the 2- 3 spectra
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Con_spectra_data_97.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "Concentration"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 2.00 and 3.00 and remove concentration 1
data_long <- data_long %>%
  filter(PPM >= 2.00 & PPM <= 3.00 & Concentration != "1")

# Adjust the vertical spacing by adding an offset based on Concentration
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    Concentration == "8" ~ 8 * 80000,
    Concentration == "4" ~ 4 * 80000,
    Concentration == "2" ~ 2 * 90000,
    Concentration == "0" ~ 0 * 400,
    TRUE ~ as.numeric(Concentration) * 100000  # Default case
  ))

# Define the order of Concentrations for the legend
concentration_order <- c("8", "4", "2", "0")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(Concentration, levels = concentration_order))) +
  geom_line(size = 0.5) +
  scale_x_reverse(breaks = seq(2, 3, by = 00.25), limits = c(3, 2)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Adjust y-axis limits
  labs(title = "The Effect of Concentration on the Spectra TE 97ms",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "Concentration (mM)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot_2_3_PPM.png")





################################
#spectra 6 all con 
#0 - 288


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Con_spectra_data_6.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "Concentration"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 1.00 and 5.00 and remove concentration 1
data_long <- data_long %>%
  filter(PPM >= 1.00 & PPM <= 5.00 & Concentration != "1")

# Adjust the vertical spacing by adding an offset based on Concentration
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    Concentration == "8" ~ 8 * 80000,
    Concentration == "4" ~ 4 * 80000,
    Concentration == "2" ~ 2 * 90000,
    Concentration == "0" ~ 0 * 400,
    TRUE ~ as.numeric(Concentration) * 100000  # Default case
  ))

# Define the order of Concentrations for the legend
concentration_order <- c("8", "4", "2", "0")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(Concentration, levels = concentration_order))) +
  geom_line(size = 0.5) +
  scale_x_reverse(breaks = seq(1, 5, by = 0.25), limits = c(5, 1)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Adjust y-axis limits
  labs(title = "The Effect of Concentration on the Spectra TE6ms ",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "Concentration (mM)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")



# from the 2- 3 spectra
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("Con_spectra_data_6.csv")

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "Concentration"),
               names_pattern = "(PPM|mag)_(\\d+)")

# Filter data to include only PPM values between 2.00 and 3.00 and remove concentration 1
data_long <- data_long %>%
  filter(PPM >= 2.00 & PPM <= 3.00 & Concentration != "1")

# Adjust the vertical spacing by adding an offset based on Concentration
data_long <- data_long %>%
  mutate(mag = mag + case_when(
    Concentration == "8" ~ 8 * 80000,
    Concentration == "4" ~ 4 * 80000,
    Concentration == "2" ~ 2 * 90000,
    Concentration == "0" ~ 0 * 400,
    TRUE ~ as.numeric(Concentration) * 100000  # Default case
  ))

# Define the order of Concentrations for the legend
concentration_order <- c("8", "4", "2", "0")

# Plot the spectra with separation
ggplot(data_long, aes(x = PPM, y = mag, color = factor(Concentration, levels = concentration_order))) +
  geom_line(size = 0.5) +
  scale_x_reverse(breaks = seq(2, 3, by = 00.25), limits = c(3, 2)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Adjust y-axis limits
  labs(title = "The Effect of Concentration on the Spectra TE 6ms",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "Concentration (mM)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))

# Save the plot
ggsave("spectra_plot_2_3_PPM.png")






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



