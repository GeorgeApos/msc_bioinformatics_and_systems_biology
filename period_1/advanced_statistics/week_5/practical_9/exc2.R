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
daphnia_data <- read.table("DaphniaData.txt", header = TRUE)
daphnia_data$temperature <- as.factor(daphnia_data$temperature)
daphnia_data$quality <- as.factor(daphnia_data$quality)

head(daphnia_data)

# Interaction plot (raw velocity)
with(daphnia_data, interaction.plot(x.factor = quality,
                                    trace.factor = temperature,
                                    response = velocity,
                                    fun = mean,
                                    type = "b",
                                    legend = TRUE,
                                    xlab = "water quality",
                                    ylab = "daphnia mean velocity",
                                    pch = c(1, 2)))

# Full factorial ANOVA on raw velocity
daphnia.lm1 <- lm(velocity ~ temperature * quality, data = daphnia_data)
summary(daphnia.lm1)
Anova_table(daphnia.lm1)

  # Diagnostics (raw model)
plot(daphnia.lm1)

# Log-transform velocity
daphnia_data$velocity <- log(daphnia_data$velocity)

# Interaction plot (log-transformed velocity)
with(daphnia_data, interaction.plot(x.factor = quality,
                                    trace.factor = temperature,
                                    response = velocity,
                                    fun = mean,
                                    type = "b",
                                    legend = TRUE,
                                    xlab = "water quality",
                                    ylab = "ln(velocity)",
                                    pch = c(1, 2)))

# ANOVA on log-transformed velocity
daphnia.lm2 <- lm(velocity ~ temperature * quality, data = daphnia_data)
summary(daphnia.lm2)
Anova_table(daphnia.lm2)

# Diagnostics (log model)
plot(daphnia.lm2)

# Pairwise comparisons (Fisher’s LSD, no adjustment)
emm1a <- emmeans(daphnia.lm2, specs = "quality")
pairs(emm1a, adjust = "none")

# Means for all treatments
emm1 <- emmeans(daphnia.lm2, specs = c("quality", "temperature"))
emm1

# Back-transform confidence intervals
semm1 <- summary(emm1)
exp(semm1[, c("emmean", "SE", "lower.CL", "upper.CL")])
  
