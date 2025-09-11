# Load data
organic_farming <- read.table("OrganicFarming.txt", 
                              header = TRUE,
                              stringsAsFactors = TRUE,
                              sep = "\t",
                              na.strings = "NA",
                              dec = ".",
                              strip.white = TRUE
)

# Quick preview
head(organic_farming)

organic_farming$Certification <- as.factor(organic_farming$Certification)

with(organic_farming,
     Barplot(Type,
             by = Certification,
             style = "divided",
             legend.pos = "above",
             xlab = "Type of produce",
             ylab = "Percent",
             scale = "percent",
             label.bars = TRUE
     )
)

# Summarize data in a table
organic_farming_table <- xtabs(formula = ~ Type + Certification, data = organic_farming)

# Print table
cat("\nFrequency table:\n")
print(organic_farming_table)

# Chi-square test
organic_farming_test <- chisq.test(organic_farming_table, correct = FALSE)
print(organic_farming_test)

# Expected counts
cat("\nExpected counts:\n")
print(organic_farming_test$expected)