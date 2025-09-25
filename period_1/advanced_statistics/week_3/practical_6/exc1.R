# Get overview of all packages on your computer
pkgs <- rownames(installed.packages())

# If ggplot2 package is missing, install it
if (!"ggplot2" %in% pkgs) install.packages("ggplot2")

library(ggplot2)

# Read data into an object named car_data
# N.B.: Reading in a header (= column names) can be TRUE or FALSE.
#       In our data files, we always provide column names, so we specify header = TRUE

car_data <- read.table("Car sales.txt", header = TRUE)

# Complete the code of the following command such that 5 data lines are shown

head(x=car_data, n=5)

# Visualise the pairwise relations between the variables

pairs(car_data[,c(6,4,5)])

# estimate Model 1

Model1 <- lm(formula = ln_price ~ engine_s, data = car_data)
summary(Model1)

# estimate Model 2

Model2 <- lm(formula = ln_price ~ horsepow, data = car_data)
summary(Model2)

# estimate Model 3

Model3 <- lm(formula = ln_price ~ engine_s + horsepow, data = car_data)
summary(Model3)

# create a bar chart showing the size of each horsepower bin
# first make sure we get the right order of our bars
# next create bar chart, with some options to make it prettier

car_data$pk_binned<-factor(car_data$pk_binned, 
                           levels=c("(54.6,81.3]", "(81.3,108]", "(108,134]","(134,160]", "(160,187]", "(187,213]",
                                    "(213,239]", "(239,266]","(266,292]", "(292,318]", "(345,371]", "(424,450]"))

barplot(table(car_data$pk_binned), las=2, cex.names=0.7, xlab="horse power", ylab="nr of models")

# create a scatter plot of ln(price) vs engine size, 
# with groups coloured according to horsepower
# and linear regression lines per group

ggplot(data = car_data, aes(y=ln_price,x=engine_s,
                            color=pk_binned)) + 
  geom_point()+ 
  stat_smooth(method="lm",se=FALSE) 

