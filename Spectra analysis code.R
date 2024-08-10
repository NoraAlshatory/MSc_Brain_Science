
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
  scale_x_reverse(breaks = seq(1.25, 4.3, by = 0.25), limits = c(4.3, 1.25)) +  # Set specific breaks for the x-axis and limits
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
  filter(PPM >= 1.25 & PPM <= 4.3)  # Adjusted filter for PPM values between 1.25 and 4.3

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
  scale_x_reverse(breaks = seq(1.25, 4.3, by = 0.25), limits = c(4.3, 1.25)) +  # Set specific breaks for the x-axis and limits
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  labs(title = "Control Spectra at Different TEs (STEAM)",
       x = "Chemical Shift (PPM)",
       y = "Magnitude",
       color = "TE (ms)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 19, face = "bold"),
        legend.text = element_text(size = 19, face = "bold"))

# Save the plot
ggsave("spectra_plot.png")







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
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.9) +
  scale_x_reverse(limits = c(4.3, 1.25), breaks = seq(1.25, 4.3, by = 0.25)) +  # Adjust x-axis limits
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  geom_vline(xintercept = 2.25, linetype = "longdash", color = "red", size = 0.5) +  # Add a vertical line at 2.25 PPM
  geom_vline(xintercept = 2.72, linetype = "longdash", color = "blue", size = 0.5) +  # Add a vertical line at 2.72 PPM
  labs(title = " (C)  Phantom 2 Separated Spectra at Different TEs (PRESS)",
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
ggsave("spectra_plot_peaks_2.25_2.72.png")







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

# Filter data to include only PPM values between 1.25 and 4.3
data_long <- data_long %>%
  filter(PPM >= 1.25 & PPM <= 4.3)

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
te_order <- c("97", "78", "35", "6")

# Plot the spectra with separation and add dashed lines at 2.25 PPM and 2.72 PPM
ggplot(data_long, aes(x = PPM, y = mag, color = factor(TE, levels = te_order))) +
  geom_line(size = 0.9) +
  scale_x_reverse(breaks = seq(1.25, 4.3, by = 0.25), limits = c(4.3, 1.25)) +  # Set specific breaks for the x-axis
  scale_y_continuous(limits = c(min(data_long$mag), max(data_long$mag) + 100000)) +  # Increase the y-axis limits
  geom_vline(xintercept = 2.25, linetype = "longdash", color = "red", size = 0.5) +  # Add a vertical line at 2.25 PPM
  geom_vline(xintercept = 2.72, linetype = "longdash", color = "blue", size = 0.5) +  # Add a vertical line at 2.72 PPM
  labs(title = " (D)  Phantom 2 Separated Spectra at Different TEs (STEAM)",
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
ggsave("spectra_plot_with_peaks.png")






