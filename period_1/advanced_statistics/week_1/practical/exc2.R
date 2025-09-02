# Read data into an object named Yield_data. 
Yield_data <- read.table(file = "Corn_yield.txt", header = TRUE)

# View data
View(Yield_data)

# Create a boxplot for yield of each variety
with(data = Yield_data, boxplot(yield ~ Variety))

# Create a strip chart for yield for each variety.
stripchart(yield ~ Variety, vertical=TRUE, method="jitter", ylab="yield", data=Yield_data)

# Perform an independent samples t-test and assume equal variances for both groups
t.test(yield ~ Variety, alternative = 'two.sided', conf.level = 0.95, var.equal = TRUE,  data = Yield_data)

# Perform Levene's test  
leveneTest(yield ~ Variety, data = Yield_data, center = "mean")