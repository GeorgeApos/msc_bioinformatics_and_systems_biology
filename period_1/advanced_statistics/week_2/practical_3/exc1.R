# Load data
sports_car <- read.table("sports cars full.txt",
                         header = TRUE,
                         stringsAsFactors = TRUE,
                         sep = "\t",
                         na.strings = "NA",
                         dec = ".",
                         strip.white = TRUE
)

# Quick preview of first few data lines
head(sports_car)

# Summarize dataset into a table
sports_car_table <- with(sports_car, table(outcome))

# Print counts
print(sports_car_table)

# Print percentages
print(100 * sports_car_table / sum(sports_car_table))

# Exact binomial test (one-sided)
binom.test(x = 90, n = 150, p = 0.7, alternative = "less", conf.level = 0.95)

# For confidence interval: Exact binomial test (two-sided)
binom.test(x = 90, n = 150, p = 0.7, alternative = "two.sided", conf.level = 0.95)

### using the table we made before
binom.test(sports_car_table, p=0.3, alternative = "greater", conf.level = 0.95)

# The one-tailed p-value using the normal approximation without continuity correction
prop.test(x = 30, n = 50, p = 0.7, alternative = "less", correct = FALSE)

