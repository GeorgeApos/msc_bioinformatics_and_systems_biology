###############################################################
## Practical 11 – Exercise 3: Spacing & Nitrogen (Two-way ANOVA + Blocks)
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
spacing_data <- read.table(file = "Spacing and nitrogen.txt", header = TRUE)
spacing_data$block <- as.factor(spacing_data$block)
spacing_data$Nitrogen <- as.factor(spacing_data$Nitrogen)
spacing_data$Spacing <- as.factor(spacing_data$Spacing)

# --- Interaction plot ignoring blocks ---
with(spacing_data, interaction.plot(
  x.factor = Nitrogen,
  trace.factor = Spacing,
  response = Yield,
  type = "b",
  pch = c(1, 16, 24)
))

# --- Fit full model ---
spacing.lm1 <- lm(Yield ~ Nitrogen + Spacing + Nitrogen:Spacing + block, data = spacing_data)
summary(spacing.lm1)
Anova_table(spacing.lm1)

# --- Profile plots (emmeans) ---
emmip(spacing.lm1, block ~ Spacing)
emmip(spacing.lm1, block ~ Nitrogen)

# --- LSD pairwise comparisons ---
pairs(emmeans(spacing.lm1, specs = "Spacing"), adjust = "none")
pairs(emmeans(spacing.lm1, specs = "Nitrogen"), adjust = "none")
