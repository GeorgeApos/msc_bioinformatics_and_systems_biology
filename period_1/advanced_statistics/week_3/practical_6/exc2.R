## Read data into an object named salary_data

salary_data <- read.table(file = "Programmer Salaries.txt", header = TRUE)

# complete the code for head() to check out the data
# head() will provide the first n rows
# summary() gives some summary info per variable

head(x = salary_data, n=10)
summary(object = salary_data)

# Create a histogram of salary

with(salary_data, hist(x = Salary))
  
# Create matrix plot for all variables, colour = blue

pairs(salary_data, col = "blue")

# Create a linear model; complete the code

salary.mlr <- lm(formula = Salary ~ Numemp + Margin + IPCost, data = salary_data)
summary(salary.mlr)

# Inspect the fit of the model
# Plot observed salaries against the predictions of salaries based on the model, colour purple

predicted <- fitted.values(salary.mlr)
plot(y = salary_data$Salary, x = predicted, col = "purple")

