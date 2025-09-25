# read the data into the object fish_data
fish_data <- read.table(file = "Fish quality vs delay.txt", header = TRUE)
# create the scatter plot as in the previous practical
with(fish_data, plot(fish_quality ~ delay, 
                     pch = 16, 
                     xlab = "Delay in ice storage/[hours]",
                     ylab = "Fish quality/[ten-point scale]"))

# create linear model; complete the code. Use lm().
fishquality.lm1 <- lm(formula = fish_quality ~ delay, data = fish_data)
summary(fishquality.lm1)

fishquality.lm2 <- lm(formula = fish_quality ~ delay + I(delay^2), data = fish_data)
summary(fishquality.lm2)

## Make a residual diagram, where the residuals are plotted versus the delay in
## ice storage in hours (x)
plot(x = fish_data$delay,
     y = fishquality.lm1$resid,
     xlab = "Delay in ice storage/[hours]",
     ylab = "Residual",
     col = "red")
# add horizontal line at 0
abline(h=0)

## data scatterplot
with(fish_data, plot(fish_quality ~ delay, 
                     pch = 16, 
                     xlab = "Delay in ice storage/[hours]",
                     ylab = "Fish quality/[ten-point scale]"))


# add the model lm1 to the plot
abline(fishquality.lm1, col = "red")

# create variable delay2 within the fish_data object (using $)
# that is delay squared; complete the code; for x-squared 
# you can use x^2 or x*x.
fish_data$delay2 <- fish_data$delay^2

# compare lm1 and lm2
anova(fishquality.lm1, fishquality.lm2)

# create again the scatter plot of quality vs delay
# abline only creates straight lines, so we use the 'lines' function
# first create predicted y-values for both models, for every x-value
pred1 <- fitted.values(fishquality.lm1)
pred2 <- fitted.values(fishquality.lm2)
# scatter-plot
with(fish_data, plot(fish_quality ~ delay, 
                     pch = 16,                             # filled cricles
                     xlab = "Delay in ice storage/[hours]",
                     ylab = "Fish quality/[ten-point scale]"))
# line connecting predicted values for model lm1
# note that the x-values are fish_data$delay
lines(pred1 ~ fish_data$delay, lty = 1, col = "red")
# line connecting predicted values for model lm2
lines(pred2 ~ fish_data$delay, lty = 2, col = "blue")

# this is rather crude; with more effort we can create a
# curved line, but we try to keep the code as simple as possible

# same as before, but now for lm2


plot(x = fish_data$delay,
     y = fishquality.lm2$resid,
     xlab = "Delay in ice storage/[hours]",
     ylab = "Residual",
     col = "red")
# add horizontal line at 0
abline(h=0)