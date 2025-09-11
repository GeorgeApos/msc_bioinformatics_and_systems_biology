# Enter data
steroids <- matrix(c(103, 8440,
                     52, 4289,
                     65, 6428),
                   nrow = 3,
                   ncol = 2,
                   byrow = TRUE)

dimnames(steroids) <- list("Division"=c("I", "II", "III"), "Admitted"=c("Yes", "No"))

# Counts
steroids

# Row Percentages
rowPercents(steroids)

# Chi-square test
chisq.test(steroids, correct = FALSE)
