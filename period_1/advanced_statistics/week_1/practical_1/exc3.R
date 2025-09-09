# Read the data in R
coal_data <- read.table("coal_mining_paired.txt", header = TRUE)

# View the dataset
View(coal_data)

# Perform a paired t-test for the difference in pH 
with(coal_data,
     t.test(before, after, alternative='two.sided', conf.level=.95, paired=TRUE))   

# Create a new variable
coal_data$d <- with(coal_data, before - after)

# Perform a one-sample t-test on d
with(coal_data, t.test(d, alternative='two.sided', conf.level = 0.95))

# Create QQ-plot to check the normality assumption of d
with(coal_data, qqPlot(d))

normalityTest( ~d, test = "shapiro.test", data = coal_data)
