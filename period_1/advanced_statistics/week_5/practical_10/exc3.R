## Load paired data
mining_paired_data <- read.table(file = "coal_mining_paired.txt", header = TRUE)
mining_paired_data$plot <- as.factor(mining_paired_data$plot)

## Inspect
head(mining_paired_data)

## Paired t-test
with(mining_paired_data, 
     t.test(x = before, y = after, paired = TRUE))

## Reshape to block design format
mining_block_data <- melt(mining_paired_data, 
                          id.vars = c("plot"),
                          variable.name = "Moment",
                          value.name = "pH")

## Inspect reformatted data
View(mining_block_data)

## RCBD model
mining_RCBD.lm <- lm(pH ~ plot + Moment, data = mining_block_data)
summary(mining_RCBD.lm)
Anova_table(mining_RCBD.lm)

## Means per Moment
summaryBy(pH ~ Moment, data = mining_block_data, FUN = mean)
