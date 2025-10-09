# ---------------------------------------------------------
# Exercise 2: Genotype and dbh
# ---------------------------------------------------------

library(car)
library(ggplot2)

# Read and prepare data
dbh_data <- read.table(file = "Genotype_dbh.txt", header = TRUE)
dbh_data$genotype <- as.factor(dbh_data$genotype)
dbh_data$thickness <- as.factor(dbh_data$thickness)
levels(dbh_data$thickness) <- c("thin", "medium", "thick")
levels(dbh_data$genotype) <- c("AB", "Ab", "aB", "ab")

# ---------------------------------------------------------
# (a) Chi-square test for expected genotype proportions
# ---------------------------------------------------------
props <- c(0.5, 0.25, 0.17, 0.08)
observed <- table(dbh_data$genotype)
chisq.test(observed, p = props)
chisq.test(observed, p = props)$expected

# ---------------------------------------------------------
# (b) Test for independence between thickness and genotype
# ---------------------------------------------------------
with(dbh_data, chisq.test(genotype, thickness))
with(dbh_data, chisq.test(genotype, thickness))$observed
with(dbh_data, chisq.test(genotype, thickness))$expected

# ---------------------------------------------------------
# (c) Graphical exploration
# ---------------------------------------------------------
with(dbh_data, plot(dbh ~ thickness))
with(dbh_data, plot(dbh ~ genotype))

ggplot(dbh_data, aes(x = genotype, y = dbh)) +
  geom_dotplot(aes(fill = thickness, color = thickness),
               binaxis = 'y', stackdir = 'center', binwidth = .3)

ggplot(dbh_data, aes(x = genotype, y = dbh)) +
  geom_boxplot(aes(color = thickness))

# ---------------------------------------------------------
# (d) One-way ANOVA for dbh differences among genotypes
# ---------------------------------------------------------
Anova_table <- function(model) {
  At2 <- Anova(model, Type = 2)
  At1 <- anova(model)
  SStotal <- sum(At1[, 2])
  dftotal <- sum(At1[, 1])
  Total <- as.data.frame(cbind(SStotal, dftotal, "", ""))
  rownames(Total) <- "Total"
  names(Total) <- names(At2)
  Avtable <- rbind(At2, Total)
  for (c in 1:4) { Avtable[, c] <- as.numeric(Avtable[, c]) }
  rows <- dim(Avtable)[1]
  MS <- as.numeric(Avtable[, 1]) / as.numeric(Avtable[, 2])
  Avtable <- add_column(Avtable, MS, .after = 2)
  return(Avtable)
}

dbh.oneway <- lm(dbh ~ genotype, data = dbh_data)
summary(dbh.oneway)
Anova_table(dbh.oneway)

# ---------------------------------------------------------
# (e) Test for equal variances and residuals
# ---------------------------------------------------------
leveneTest(dbh.oneway, center = mean)
plot(dbh.oneway)  # diagnostic plots
