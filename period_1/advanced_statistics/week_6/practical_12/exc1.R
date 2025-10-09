# ---------------------------------------------------------
# Exercise 1: Type and amount of soap
# ---------------------------------------------------------

# Load required packages
pkgs <- rownames(installed.packages())
if (!"car" %in% pkgs) install.packages("car")
if (!"ggplot2" %in% pkgs) install.packages("ggplot2")
if (!"fastDummies" %in% pkgs) install.packages("fastDummies")

library(car)
library(ggplot2)
library(fastDummies)
library(tibble)

# Define Anova_table function
Anova_table <- function(model) {
  At2 <- Anova(model, Type = 2)
  At1 <- anova(model)
  SStotal <- sum(At1[, 2])
  dftotal <- sum(At1[, 1])
  Total <- as.data.frame(cbind(SStotal, dftotal, "", ""))
  rownames(Total) <- "Total"
  names(Total) <- names(At2)
  Avtable <- rbind(At2, Total)
  for (c in 1:4) { Avtable[, c] <- as.numeric(Avtable[, c]) }
  rows <- dim(Avtable)[1]
  MS <- as.numeric(Avtable[, 1]) / as.numeric(Avtable[, 2])
  Avtable <- add_column(Avtable, MS, .after = 2)
  return(Avtable)
}

# ---------------------------------------------------------
# Data and model fitting
# ---------------------------------------------------------

# Read data
soap_data <- read.table(file = "Laundry.txt", header = TRUE)
soap_data$Type <- as.factor(soap_data$Type)

# Model 1: parallel lines (no interaction)
soap.lm1 <- lm(soap_height ~ soap_amount + Type, data = soap_data)
summary(soap.lm1)
Anova_table(soap.lm1)

# Plot fitted values
ggplot(soap_data, aes(x = soap_amount, y = soap_height, color = Type)) +
  geom_point() +
  geom_line(aes(y = soap.lm1$fitted.values), size = 1)

# Model 2: interaction model (different slopes)
soap.lm2 <- lm(soap_height ~ soap_amount * Type, data = soap_data)
summary(soap.lm2)
Anova_table(soap.lm2)

# Compare models
anova(soap.lm1, soap.lm2)
