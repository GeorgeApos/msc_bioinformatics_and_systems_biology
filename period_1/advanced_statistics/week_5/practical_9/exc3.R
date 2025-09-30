# Load libraries
library(car)
library(tibble)
library(ggplot2)
library(emmeans)

# Function for Type II ANOVA table
Anova_table <- function(model) {
  At2 <- Anova(model, Type = 2)         # Type II SS
  At1 <- anova(model)                   # Type I SS
  SStotal <- sum(At1[, 2])
  dftotal <- sum(At1[, 1])
  Total <- as.data.frame(cbind(SStotal, dftotal, "", ""))
  rownames(Total) <- "Total"
  names(Total) <- names(At2)
  Avtable <- rbind(At2, Total)
  
  for (c in 1:4) {
    Avtable[, c] <- as.numeric(Avtable[, c])
  }
  rows <- dim(Avtable)[1]
  MS <- as.numeric(Avtable[, 1]) / as.numeric(Avtable[, 2])
  Avtable <- add_column(Avtable, MS, .after = 2)
  return(Avtable)
}

# Load data
orange_data <- read.table("OrangeTreeTrunks.txt", header = TRUE)
orange_data$pH <- as.factor(orange_data$pH)
orange_data$Calcium <- as.factor(orange_data$Calcium)

# Full factorial ANOVA
orange.lm1 <- lm(Trunkdiameter ~ pH * Calcium, data = orange_data)
summary(orange.lm1)
Anova_table(orange.lm1)

# Profile plots
emmip(orange.lm1, pH ~ Calcium)
emmip(orange.lm1, Calcium ~ pH)
