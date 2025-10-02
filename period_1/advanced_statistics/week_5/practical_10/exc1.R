## Install packages if missing
pkgs <- rownames(installed.packages())

if (!"car" %in% pkgs) install.packages("car")
if (!"doBy" %in% pkgs) install.packages("doBy")
if (!"dplyr" %in% pkgs) install.packages("dplyr")
if (!"emmeans" %in% pkgs) install.packages("emmeans")
if (!"ggplot2" %in% pkgs) install.packages("ggplot2")
if (!"reshape2" %in% pkgs) install.packages("reshape2")
if (!"tibble" %in% pkgs) install.packages("tibble")

## Load libraries
library(car)
library(doBy
library(dplyr)
library(emmeans)
library(ggplot2)
library(reshape2)
library(tibble)

## Custom ANOVA function (Type II SS + total row)
Anova_table <- function(model) {
  At2 <- Anova(model, Type = 2)     # Type II ANOVA
  At1 <- anova(model)               # Type I ANOVA
  SStotal <- sum(At1[, 2])
  dftotal <- sum(At1[, 1])
  
  Total <- as.data.frame(cbind(SStotal, dftotal, "", ""))
  rownames(Total) <- "Total"
  names(Total) <- names(At2)
  
  Avtable <- rbind(At2, Total)
  
  for (c in 1:4) {
    Avtable[, c] <- as.numeric(Avtable[, c])
  }
  
  rows <- dim(Avtable)[1]
  MS <- as.numeric(Avtable[, 1]) / as.numeric(Avtable[, 2])
  Avtable <- add_column(Avtable, MS, .after = 2)
  
  return(Avtable)
}

## Load data
music_data <- read.table(file = "ex15-6.txt", header = TRUE)
music_data$TypeMusic <- as.factor(music_data$TypeMusic)
music_data$Subject   <- as.factor(music_data$Subject)

## Inspect
str(music_data)

## Descriptive statistics
music_data %>% 
  group_by(TypeMusic) %>% 
  summarise(
    mean = sprintf("%0.2f", mean(Efficiency)),
    sd   = sprintf("%0.2f", sd(Efficiency)),
    count = n()
  )

## Boxplot
ggplot(data = music_data, aes(x = TypeMusic, y = Efficiency)) +
  geom_boxplot(fill = "darkgrey") +
  theme_bw()

## Fit RCBD model (additive, no interaction)
musicRCBD.lm1 <- lm(Efficiency ~ TypeMusic + Subject, data = music_data)
summary(musicRCBD.lm1)
Anova_table(musicRCBD.lm1)

## Interaction plot
with(music_data, interaction.plot(x.factor = TypeMusic, 
                                  trace.factor = Subject, 
                                  response = Efficiency))

## Full model with interaction (overfitted!)
musicRCBD.lm0 <- lm(Efficiency ~ TypeMusic * Subject, data = music_data)
summary(musicRCBD.lm0)

# Uncomment below to see warnings
# Anova_table(musicRCBD.lm0)
# anova(musicRCBD.lm0)

## Estimated Marginal Means
emm <- emmeans(musicRCBD.lm1, specs = "TypeMusic")
emm

## Tukey post-hoc test
pairs(emm, adjust = "Tukey")
