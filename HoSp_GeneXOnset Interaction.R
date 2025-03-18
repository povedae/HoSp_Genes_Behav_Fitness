#Install Libraries
library(dplyr)
library(lubridate)
library(ggplot2)
library(readr)
library(emmeans)
library(multcomp)

#Load data
genes <- read.csv("hs.csv")
onset <- read.csv("onset.csv")

# Merge the datasets
merged_data <- inner_join(genes, onset, by = c("ID", "site"))


# Filter the data for ZT == 1
ZT1_data <- merged_data %>%
  filter(ZT == 1 & Gene == 'BMAL')

# Define the shapes for each site
site_shapes <- c("farm" = 15, "greens" = 16, "campus" = 17)  # square, circle, triangle

# Plot activity onset for ZT == 1
ggplot(ZT1_data, aes(x = Expression, y = avg_onset_rel, color = site, shape = site)) +
  geom_point(size = 4) +
  geom_smooth(aes(group = site), method = "lm", se = TRUE) +
  scale_colour_manual(values = c("farm" = "#5f5953", "greens" = "#1A25FF", "campus" = "#F2EE05")) +
  scale_shape_manual(values = site_shapes) +
  labs(title = "BMAL at ZT = 1",
       x = "Bmal1 Expression at ZT1",
       y = "Mean relative activity onset (mins)",
       color = "Site",
       shape = "Site") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 14) 
  ) +
  ylim(min(ZT1_data$avg_onset_rel), 25)  # Set y-axis limit

# Regression model with interaction term
regression_model <- lm(Expression ~ avg_onset_rel * site, data = ZT1_data)
summary(regression_model)


# Get R-squared and p-values for each site separately
site_list <- unique(ZT1_data$site)

# Initialize an empty list to store results
site_stats <- list()

for (site in site_list) {
  # Filter data for the specific site
  site_data <- ZT1_data %>% filter(site == !!site)
  
  # Fit the linear model
  model <- lm(Expression ~ avg_onset_rel, data = site_data)
  
  # Extract R-squared value
  r_squared <- summary(model)$r.squared
  
  # Extract p-value for the model
  p_value <- summary(model)$coefficients[2, 4]  # p-value for avg_onset_rel
  
  # Store results
  site_stats[[site]] <- list(R_squared = r_squared, P_value = p_value)
}

# Print R-squared and p-values for each site
site_stats


#Post-hoc comparison

# Compute estimated marginal means
emm <- emmeans(regression_model, pairwise ~ site, adjust = "tukey")

# Display results
emm


