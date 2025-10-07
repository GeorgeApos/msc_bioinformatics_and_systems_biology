###############################################################
## Practical 11 – Exercise 1: Scientific Writing (ANCOVA)
## Advanced Statistics (MAT20306)
###############################################################

# Install/load required packages
pkgs <- rownames(installed.packages())
if (!"car" %in% pkgs) install.packages("car")
if (!"emmeans" %in% pkgs) install.packages("emmeans")
if (!"tibble" %in% pkgs) install.packages("tibble")
if (!"ggplot2" %in% pkgs) install.packages("ggplot2")
if (!"doBy" %in% pkgs) install.packages("doBy")

library(car)
library(emmeans)
library(tibble)
library(ggplot2)
library(doBy)

# Custom ANOVA table function
Anova_table <- function(model) {
  At2 <- Anova(model, Type = 2)
  At1 <- anova(model)
  SStotal <- sum(At1[, 2])
  dftotal <- sum(At1[, 1])
  Total <- as.data.frame(cbind(SStotal, dftotal, "", ""))
  rownames(Total) <- "Total"
  names(Total) <- names(At2)
  Avtable <- rbind(At2, Total)
  for (c in 1:4) Avtable[, c] <- as.numeric(Avtable[, c])
  rows <- dim(Avtable)[1]
  MS <- as.numeric(Avtable[, 1]) / as.numeric(Avtable[, 2])
  Avtable <- add_column(Avtable, MS, .after = 2)
  return(Avtable)
}

# --- Load data ---
teaching_data <- read.table(file = "Teaching material wur.txt", header = TRUE)
teaching_data$Program <- as.factor(teaching_data$Program)

# --- Plot data ---
ggplot(data = teaching_data, aes(y = score, x = age, color = Program)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)

# --- Fit ANCOVA model ---
teaching.lm <- lm(score ~ Program + age, data = teaching_data)
summary(teaching.lm)
Anova_table(teaching.lm)

# --- Plot fitted ANCOVA model ---
ggplot(teaching_data, aes(x = age, y = score, color = Program)) +
  geom_point() +
  geom_line(aes(y = teaching.lm$fitted.values), linewidth = 1)

# --- Check mean ages ---
summaryBy(age ~ Program, data = teaching_data, FUN = mean)
mean(teaching_data$age)

# --- Pairwise comparisons (LSD) ---
emm1 <- emmeans(teaching.lm, specs = "Program", by = "age")
pairs(emm1, adjust = "none")
