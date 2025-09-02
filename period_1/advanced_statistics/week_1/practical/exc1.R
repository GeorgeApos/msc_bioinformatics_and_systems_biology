# Get overview of all packages on your computer
pkgs <- rownames(installed.packages())

# If one of our packages is missing, install it
if (!"RcmdrMisc" %in% pkgs) install.packages("RcmdrMisc")
if (!"ggplot2" %in% pkgs) install.packages("ggplot2")
if (!"car" %in% pkgs) install.packages("car")

library(RcmdrMisc)
library(ggplot2)
library(car)

# Read data into an object named Oxygen_data. 
Oxygen_data <- read.table(file = "Oxygen.txt", header = TRUE)

# View data
View(Oxygen_data)

# Visualize data using a boxplot
with(data = Oxygen_data, boxplot(Oxygen))

# Numerical summary of Oxygen
numSummary(Oxygen_data[ , "Oxygen", drop = FALSE],
           statistics=c("mean", "sd", "se(mean)", "quantiles"),
           quantiles=c(0, 0.25, 0.5, 0.75, 1))

# Confidence interval for the mean 
with(Oxygen_data,
     t.test(Oxygen, alternative = 'two.sided', mu = 0.0, conf.level = 0.95))

# Do a single sample t-test to check whether oxygen level is less than 5   
with(Oxygen_data,
     t.test(Oxygen, alternative = 'less', mu = 5, conf.level = 0.95))


