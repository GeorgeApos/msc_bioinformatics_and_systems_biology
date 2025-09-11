# Enter data
tumour_data <- matrix(c(10, 90,
                        14, 86,
                        19, 81),
                      nrow = 3,
                      ncol = 2,
                      byrow = TRUE
)

dimnames(tumour_data) <- list("Treatment" = c("Control", "Low dose", "High dose"),
                              "Tumours" = c("Yes", "No")
)

# Counts
tumour_data

# Row Percentages
rowPercents(tumour_data)

## Binomial test
# Control group
binom.test(x = 10, n = 100, p = 0.08, alternative = "greater", conf.level = 0.95)
# Only for confidence interval
binom.test(x = 10, n = 100, p = 0.08, alternative = "two.sided", conf.level = 0.95)

# Low dose group
binom.test(x = 14, n = 100, p = 0.08, alternative = "greater", conf.level = 0.95)
# Only for confidence interval
binom.test(x = 14, n = 100, p = 0.08, alternative = "two.sided", conf.level = 0.95)

# High dose group
binom.test(x = 19, n = 100, p = 0.08, alternative = "greater", conf.level = 0.95)
# Only for confidence interval
binom.test(x = 19, n = 100, p = 0.08, alternative = "two.sided", conf.level = 0.95)

# Chi-square test
tumour_chisq <- chisq.test(tumour_data, correct = FALSE)
print(tumour_chisq)

# Expected counts
tumour_chisq$expected