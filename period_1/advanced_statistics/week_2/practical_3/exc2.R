# Enter data
chickenegg <- matrix(c(44, 16,
                       23, 17),
                     nrow = 2,
                     ncol = 2,
                     byrow = TRUE)

# Give names to rows and columns
dimnames(chickenegg) <- list("Supplement" = c("Oyster", "Lime"), "Result" = c("Pass", "Fail"))

# Counts
chickenegg

# Fisher's test
fisher.test(chickenegg)

# Chi-square test
chisq.test(chickenegg, correct = FALSE)

# Fisher's test
fisher.test(chickenegg, alternative = "greater")

power.prop.test(p1 = 0.4, p2 = 0.25, power = 0.8, alternative = "one.sided") 