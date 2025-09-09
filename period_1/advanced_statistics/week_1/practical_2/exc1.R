# Get overview of all packages on your computer
pkgs <- rownames(installed.packages())

# If one of our packages is missing, install it
if (!"car" %in% pkgs) install.packages("car")
if (!"doBy" %in% pkgs) install.packages("doBy")

library(car)
library(doBy)

# Compute sample size
power.t.test(delta = 5, sd = 4.1, power = 0.8, type = "two.sample")

# Compute power
power.t.test(n = 10, delta = 5, sd = 4.1, type = "two.sample")
