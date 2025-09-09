# Read in the data
charity_data <- read.table(file = "Charity_money.txt", header = TRUE)

# View the data
View(charity_data)

# Look at a summary of the data
summary(charity_data)

# Summary statistics for donation in both groups seperately
summary_by(charity_data, donation ~ LivingArea)

# Boxplot by groups
boxplot(donation ~ LivingArea, data = charity_data)

# QQ-plots for both types of area
with(charity_data, qqPlot(donation ~ LivingArea))