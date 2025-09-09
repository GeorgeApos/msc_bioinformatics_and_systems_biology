# Read in the data 
iodine_bread_data <- read.table(file="Iodine_in_Breadproducts.txt", header = TRUE)

# View the first lines of the data set
head(iodine_bread_data)

# View the data set in a separate tab
View(iodine_bread_data)

# Calculate the differences between each pair for 2011 and 2016
iodine_bread_data$diff <- iodine_bread_data$iodine2016 - iodine_bread_data$iodine2011

# Look at a summary of the data
summary(iodine_bread_data)

# Construct a QQ-plot
qqPlot(iodine_bread_data$diff)

# Shapiro-Wilk test
shapiro.test(iodine_bread_data$diff)

# Wilcoxon signed rank test for paired data
wilcox.test(iodine_bread_data$iodine2011,
            iodine_bread_data$iodine2016,
            paired=TRUE,
            alternative="greater",
            conf.level=0.95)

# Paired t-test
t.test(iodine_bread_data$iodine2011,
       iodine_bread_data$iodine2016,
       paired = TRUE,
       alternative = "greater")

# 