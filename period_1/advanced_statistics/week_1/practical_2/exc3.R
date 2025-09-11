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

# Carry out Levene's test
leveneTest(donation ~ LivingArea, data = charity_data, center = mean)

# Two-sample t-test to compare means
t.test(donation ~ LivingArea, data = charity_data, var.equal=TRUE)

# Run the Wilcoxon rank sum test
# default value for alternative = "two.sided"
wilcox.test(donation ~ LivingArea, data = charity_data)

## get the parameter estimates
Tapply(donation ~ LivingArea, median, na.action=na.omit, data=charity_data) # medians by group