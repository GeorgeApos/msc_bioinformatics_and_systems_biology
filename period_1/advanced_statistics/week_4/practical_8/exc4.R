## ------------------------------------------------------------
## Exercise 4: Phosphorus in tree leaves (apple varieties)
## ------------------------------------------------------------

## 0) Packages
# install.packages("doBy")
# install.packages("emmeans")
# install.packages("car")

library(doBy)
library(emmeans)
library(car)

## If your helper file with Anova_table() exists, source it once:
## source("Function Anova_table type 2.R")

## 1) Read in the data
phosphorus_data <- read.table(file = "P in tree leaves.txt", header = TRUE)

## Quick look
cat("\n--- Head of data ---\n")
print(head(phosphorus_data))
str(phosphorus_data)

## 2) Descriptive statistics: means per variety
cat("\n--- Sample group means per variety ---\n")
print(summaryBy(phosp ~ variety, data = phosphorus_data, FUN = mean))

## If you want full summary stats (mean, sd, se, n):
summaryBy(phosp ~ variety, data = phosphorus_data,
          FUN = function(x) c(m = mean(x), s = sd(x), se = sd(x)/sqrt(length(x)), n=length(x)))

## 3) Fit linear model (ANOVA)
phosphorus.model <- lm(formula = phosp ~ variety, data = phosphorus_data)

cat("\n--- Model summary (cornerstone coding) ---\n")
print(summary(phosphorus.model))

cat("\n--- Type II ANOVA table ---\n")
print(Anova_table(phosphorus.model))

## 4) Post-hoc LSD (Fisher’s protected LSD, adjust = "none")
emm <- emmeans(phosphorus.model, specs = "variety")

cat("\n--- Estimated marginal means (EMMs) ---\n")
print(emm)

cat("\n--- Pairwise differences (LSD, no adjustment) ---\n")
print(pairs(emm, adjust = "none"))

## Plot EMMs with 95% CI
plot(emm, horizontal = FALSE,
     xlab = "Phosphorus content", ylab = "Variety",
     main = "Estimated Marginal Means (95% CI)")

## 5) Standard error of differences (SED) check
cat("\n--- Note: SED values for differences can be read from emmeans output.\n")
cat("Formula: sed = sqrt( (2 * MSE) / n ), for equal sample sizes.\n")
cat("Example: Gala vs Jonagold.\n")

## 6) Non-parametric alternative: Rank-ANOVA
phosphorus_data$phosp_rank <- rank(phosphorus_data$phosp)
phosp_rank.lm <- lm(phosp_rank ~ variety, data = phosphorus_data)

cat("\n--- Type II ANOVA on ranks ---\n")
print(Anova_table(phosp_rank.lm))

## 7) Compare with original ANOVA
cat("\n--- Type II ANOVA on original phosp values ---\n")
print(Anova_table(phosphorus.model))

## ---- End of script ----
