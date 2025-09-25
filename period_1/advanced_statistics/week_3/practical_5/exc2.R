# Read in the data from the txt file
fishquality_data<-read.table(file="Fish quality vs delay.txt", header = TRUE)
#  check out first 5 lines of the data
head(x=fishquality_data, n=5)

# plot quality vs delay (hrs)
plot(x = fishquality_data$delay,
     y = fishquality_data$fish_quality,
     xlab = "Delay in ice storage [hours]",
     ylab = "Fish quality [ten-point scale]",
     pch = 16) # use filled circles; pch = plotting character

# Build a Simple Linear Regression model of quality vs. delay and print it; 

fishquality.lm <- lm(formula = fish_quality ~ delay, data = fishquality_data)

#  show the most important part of the results of the regression analysis
summary(object = fishquality.lm)

# Add the simple linear regression line (dashed (lty) and blue (col) )
plot(x = fishquality_data$delay,
     y = fishquality_data$fish_quality,
     xlab = "Delay in ice storage [hours]",
     ylab = "Fish quality [ten-point scale]",
     pch = 16) 
abline(reg = fishquality.lm, lty = 2, col = "blue")

# get the 0.95 confidence interval for the slope ; 
confint(object = fishquality.lm, level=0.95)

# diagnostic plots for the regression model; 
# Press Enter in the Console to see all 4 plots one by one.
plot(fishquality.lm)

# diagnostic plot : residuals vs x 
plot(fishquality.lm$residuals ~ fishquality_data$delay)
abline(0,0)
