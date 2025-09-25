## ------------------------------------------------------------
## Exercise 3: Cholesterol ~ APOA4 genotype (one-way ANOVA)
## ------------------------------------------------------------

## 0) Packages (install once if needed)
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("car")

library(ggplot2)
library(emmeans)
library(car)

## If your course helper is in a file, make sure it's sourced once per session:
## source("Function Anova_table type 2.R")   # ensures Anova_table() is available

## 1) Read data & prepare factors/labels
cholesterol_data <- read.table(file = "genetics cholesterol.txt", header = TRUE)
cholesterol_data$apoa4gen <- as.factor(cholesterol_data$apoa4gen)
levels(cholesterol_data$apoa4gen) <- c("1/1", "1/2", "2/2")

## Quick data check
cat("\n--- Data summary ---\n")
print(summary(cholesterol_data))
cat("\nCounts per genotype:\n")
print(table(cholesterol_data$apoa4gen))

## 2) Graphs
## Base plot (factor on x, response on y)
with(cholesterol_data, plot(chol ~ apoa4gen,
                            main = "Cholesterol by APOA4 Genotype",
                            xlab = "APOA4 Genotype", ylab = "Serum Cholesterol"))

## ggplot: boxplot + points
ggplot(cholesterol_data, aes(x = apoa4gen, y = chol, color = apoa4gen)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3, position = position_jitter(width = 0.08, height = 0)) +
  theme_minimal() +
  labs(title = "Serum Cholesterol by APOA4 Genotype",
       x = "APOA4 Genotype", y = "Serum Cholesterol") +
  guides(color = "none")

## 3) One-way ANOVA with lm()
cholesterol.model <- lm(formula = chol ~ apoa4gen, data = cholesterol_data)

cat("\n--- Model summary (t-tests on coefficients, cornerstone coding) ---\n")
print(summary(cholesterol.model))

cat("\n--- Type II ANOVA table ---\n")
print(Anova_table(cholesterol.model))  ## From your course helper

## 4) Estimated marginal means (group means under equal-variance model) + LSD post-hoc
emm <- emmeans(cholesterol.model, specs = "apoa4gen")
cat("\n--- Estimated marginal means (EMMs) ---\n")
print(emm)

cat("\n--- Post-hoc pairwise differences (Fisher's protected LSD: no adjustment) ---\n")
print(pairs(emm, adjust = "none"))

## 5) Plot EMMs with 95% CIs
plot(emm, horizontal = FALSE,
     comparisons = FALSE,  # set TRUE if you want arrows between pairs
     xlab = "Serum Cholesterol", ylab = "APOA4 Genotype",
     main = "Estimated Marginal Means (95% CI)")

## 6) Model diagnostics (assumptions)
cat("\n--- Diagnostic plots: Press <Enter> to cycle through ---\n")
plot(cholesterol.model)
# Inspect:
# 1) Residuals vs Fitted: look for no pattern & roughly equal spread (homoscedasticity)
# 2) Normal Q-Q: points ~ straight line (normality of errors)

## 7) Levene's test for equal variances (robust to non-normality)
cat("\n--- Levene's Test for Homogeneity of Variance ---\n")
print(leveneTest(chol ~ apoa4gen, data = cholesterol_data))
# Interpretation: p-value > 0.05 -> do not reject equal variances (assumption OK)

## 8) Rank-ANOVA (nonparametric alternative if normality is doubtful)
chol.rank.mod <- lm(formula = rank(chol) ~ apoa4gen, data = cholesterol_data)

cat("\n--- Type II ANOVA on ranks ---\n")
print(Anova_table(chol.rank.mod))

cat("\n--- Type II ANOVA on original data (for comparison) ---\n")
print(Anova_table(cholesterol.model))

## 9) Welch's ANOVA (does not assume equal variances) - FYI
cat("\n--- Welch's ANOVA (unequal variances allowed) ---\n")
print(oneway.test(chol ~ apoa4gen, data = cholesterol_data, var.equal = FALSE))

## ---- End of script ----
