# Get overview of all packages on your computer
pkgs <- rownames(installed.packages())

# If one of our packages is missing, install it
if (!"abind" %in% pkgs) install.packages("abind")
if (!"RcmdrMisc" %in% pkgs) install.packages("RcmdrMisc")

library(abind, pos = 17)
library(RcmdrMisc)

# Load data
mekong_fish <- read.table("Mekong_fish.txt",
                          header = TRUE,
                          stringsAsFactors = TRUE,
                          sep = "\t",
                          na.strings = "NA",
                          dec = ".",
                          strip.white = TRUE
)

# Quick preview
head(mekong_fish)

# Summarize dataset into a table
mekong_fish_table <- with(mekong_fish, table(species))

# Print counts
print(mekong_fish_table)

# Print percentages
print(round(100 * mekong_fish_table / sum(mekong_fish_table), 2))

# Chi-square goodness-of-fit test
# Species are in alphabetic order in R
# Thus the probabilities are not entered in the same order
# as listed in the text above, but alphabetically
mekong_fish_probs <- c(0.30, 0.50, 0.10, 0.10)
chisq.test(mekong_fish_table, p = mekong_fish_probs)
