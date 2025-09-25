# Get overview of all packages on your computer
pkgs <- rownames(installed.packages())

# If one of our packages is missing, install it
if (!"car" %in% pkgs) install.packages("car")
if (!"ggplot2" %in% pkgs) install.packages("ggplot2")
if (!"fastDummies" %in% pkgs) install.packages("fastDummies")

library(car)
library(ggplot2)
library(fastDummies)

salary_data <- read.table(file = "Programmer Salaries.txt", header = TRUE)
head(salary_data)
summary(salary_data)

# estimate the linear full model; complete the code; 
# note that we must give lm() a formula and the data-object
salaries_lm.full <- lm(formula = Salary ~ Numemp + Margin + IPCost, data = salary_data)
summary(salaries_lm.full)

## compute Pearson correlations between all variables
cor(salary_data, method="pearson")

## to get variance inflation factors, we need the package 'car'
vif(salaries_lm.full)

# estimate reduced model
salaries_lm.reduced <- lm(Salary ~ IPCost, data = salary_data)

summary(salaries_lm.full)
summary(salaries_lm.reduced)

# compare models
anova(salaries_lm.reduced, salaries_lm.full)

# critical F for rejection region
qf(p=0.05, df1=2, df2=63, lower.tail = FALSE)

