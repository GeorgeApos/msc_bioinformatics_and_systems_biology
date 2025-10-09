###############################################################
## Practical 11 – Exercise 2: Sardinella (Interaction model)
## Advanced Statistics (MAT20306)
###############################################################

# Load required packages
library(car)
library(emmeans)
library(tibble)
library(ggplot2)
library(doBy)

# Define custom ANOVA function
Anova_table <- function(model) {
  At2 <- Anova(model, Type = 2)
  At1 <- anova(model)
  SStotal <- sum(At1[, 2])
  dftotal <- sum(At1[, 1])
  Total <- as.data.frame(cbind(SStotal, dftotal, "", ""))
  rownames(Total) <- "Total"
  names(Total) <- names(At2)
  Avtable <- rbind(At2, Total)
  for (c in 1:4) Avtable[, c] <- as.numeric(Avtable[, c])
  rows <- dim(Avtable)[1]
  MS <- as.numeric(Avtable[, 1]) / as.numeric(Avtable[, 2])
  Avtable <- add_column(Avtable, MS, .after = 2)
  return(Avtable)
}

# --- Load data ---
fish_data <- read.table("sardinella.txt", sep = "\t", header = TRUE)
fish_data$environment <- as.factor(fish_data$environment)

# --- Plot raw data ---
ggplot(data = fish_data, aes(x = log_weight, y = log_O2, color = environment)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_grey() +
  theme_bw()

# --- Full model (different slopes + intercepts) ---
fish.lm1 <- lm(log_O2 ~ environment * log_weight, data = fish_data)
summary(fish.lm1)
Anova_table(fish.lm1)

# --- Reduced model (no environment effect) ---
fish.lm.rm <- lm(log_O2 ~ log_weight, data = fish_data)
summary(fish.lm.rm)
Anova_table(fish.lm.rm)

# --- Compare models (F-test) ---
anova(fish.lm.rm, fish.lm1)

# --- Plot without environment effect ---
ggplot(data = fish_data, aes(x = log_weight, y = log_O2)) +
  geom_point(size = 3, aes(color = environment)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_grey() +
  theme_bw()
