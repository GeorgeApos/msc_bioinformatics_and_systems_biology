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
tilapia_data <- read.table("Tilapia_growth.txt", header = TRUE, sep = "\t")

# Convert to factors
tilapia_data$diet <- as.factor(tilapia_data$diet)
tilapia_data$supplement <- as.factor(tilapia_data$supplement)

# Boxplot
ggplot(data = tilapia_data, aes(x = diet, y = WG_total, fill = supplement)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.5, end = 0.8) +
  theme_bw() +
  ylab("total weight gain (g)") +
  xlab("Treatment")

# Full interaction model
model1 <- lm(WG_total ~ diet * supplement, data = tilapia_data)
Anova_table(model1)
summary(model1)

# Estimated marginal means
emmeans(model1, specs = "diet", by = "supplement")

  # Additive model
model.add <- lm(WG_total ~ diet + supplement, data = tilapia_data)
emmip(emmeans(model.add, specs = "diet", by = "supplement"),
      supplement ~ diet)
