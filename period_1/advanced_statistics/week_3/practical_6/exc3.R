# Read data into an object named NO2_data

NO2_data <- read.table(file="Norway NO2.txt", header = TRUE)

# Check out the data

head(NO2_data)

# Check correlations between all variables
# To limit the output to 2 digits, we added the function round() in the second line

cor(NO2_data)
round(cor(NO2_data),2)

# Visualise the data in two sets:
# 1) y, x1, x2, x3
# 2) y, x4, x5, x6, x7
# In this case we are using column numbers in the code, rather than column names. That is just a choice

pairs(NO2_data[,1:4])
pairs(NO2_data[,c(1,5:8)])

# Estimate a linear model with all seven predictors
# Extract the information from the model that you need using summary()

NO2_full.mlr<-lm(formula = ... ~ car_per_hour + temperature + wind_speed +
                   temp_diff + wind_dir +hour + day , data = NO2_data )
summary(...)

# Create diagnostic plots for the full model

plot(NO2_full.mlr)