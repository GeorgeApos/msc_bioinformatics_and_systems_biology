## ------------------------------------------------------------
## Exercise 2: Root respiration (Douglas spruce, SO2 gas levels)
## ------------------------------------------------------------

## Load packages (install first if needed)
# install.packages("ggplot2")  # if not installed yet

library(ggplot2)

## 2.a Read the data
# Make sure RootRespiration.txt is in your working directory
respiration_data <- read.table(file = "RootRespiration.txt", header = TRUE)

# Convert SO2gas to a factor (qualitative predictor)
respiration_data$SO2gas <- as.factor(respiration_data$SO2gas)

# Inspect data
head(respiration_data)
str(respiration_data)

## 2.b Fit ANOVA model (cornerstone representation with 25 ppm as reference)
respiration.model <- lm(rootresp ~ SO2gas, data = respiration_data)

# Summary output
summary(respiration.model)

# If you have the helper function Anova_table() from your course:
# source("Function Anova_table type 2.R")   # run this first if needed
Anova_table(respiration.model)

## 2.c Get estimate of σ^2 (pooled variance) and error df
# Residual variance estimate is from Mean Sq of residuals
anova_out <- anova(respiration.model)
sigma2_hat <- anova_out["Residuals","Mean Sq"]
sigma2_hat   # estimate of σ^2
dfE <- anova_out["Residuals","Df"]
dfE           # error degrees of freedom

## 2.d Hypothesis test (overall effect of SO2 on respiration)
# Null hypothesis: μ1 = μ2 = μ3 (no SO2 effect)
# Already tested in Anova_table(respiration.model)
# Look at F statistic and p-value in the output

## 2.e Plot all datapoints (dotplot instead of boxplot)
ggplot(respiration_data, aes(x = SO2gas, y = rootresp)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2) +
  theme_minimal() +
  labs(title = "Root respiration under SO2 treatments",
       x = "SO2 concentration (ppm)",
       y = "Root respiration (ml CO2/day)")

## 2.f Check assumptions: residual diagnostic plots
plot(respiration.model)
# --> Use the first two plots:
# 1. Residuals vs Fitted (checks equal variance, linearity)
# 2. Normal Q-Q plot (checks normality)
