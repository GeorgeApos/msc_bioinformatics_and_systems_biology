# Read the data into a data-frame; first line contains column headers
MScData <- read.table(file="MSc2020.txt", header=TRUE)

# check out the first 5 lines of data
head(x=MScData, n=5)

# Visualise the relations between age, shoesize, height, and working experience (we)
pairs(x = MScData[,c("age", "height", "shoesize", "we")])

# find the Pearson correlation coefficients between these variables
cor(x = MScData[,c("age", "height", "shoesize", "we")], method="pearson")

# are these correlations different from zero? 
with(MScData, cor.test(age, we, method = "pearson"))
with(MScData, cor.test(height, shoesize, method = "pearson"))

# Colour dots by sex to illustrate how to make such graphs
pairs(MScData[,c("age", "height", "shoesize", "we")], col=MScData$gender)

## Plot shoesize vs height
with(MScData, plot(shoesize ~ height, 
                   main = "MSc students height vs shoesize", 
                   xlab="height (cm)", 
                   ylab="shoesize (European)",
                   pch=16) )  # use filled circles rather than open circles

cor(MScData$height, MScData$shoesize, method="pearson")
cor(MScData$height, MScData$shoesize, method="spearman")

# linear least squares regression shoesize by height
# Note the structure for formula =. 
# It is always: dependent ~ independent 
lm_shoesize1 <- lm(formula = shoesize ~ height, data = MScData)
summary(lm_shoesize1)

lm_shoesize1short <- lm(shoesize ~ height, data = MScData)
summary(lm_shoesize1short)

# 0.95 confidence interval
confint(object = lm_shoesize1, level = 0.95)

# 0.95 confidence interval for mean shoe size at height of 175
# First create the object newHeight to use it in the predict-function
newHeight<-data.frame(height = 175)   
predict(object = lm_shoesize1,
        newdata = newHeight,
        se.fit = TRUE,
        interval = "confidence",
        level= 0.95)

# diagnostic plots for regression model
# (1) residuals vs fitted values
plot(lm_shoesize1$resid ~ lm_shoesize1$fitted)
# (2) qqplot of residuals
qqnorm(lm_shoesize1$resid)
