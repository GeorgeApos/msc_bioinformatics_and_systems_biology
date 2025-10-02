## Load data
potato_data <- read.table(file = "Potato nitrogen blocks.txt", header = TRUE)
potato_data$Nitrogen <- as.factor(potato_data$Nitrogen)
potato_data$block    <- as.factor(potato_data$block)

## Inspect
head(potato_data, 5)

## Interaction plot (profile plot)
with(potato_data, interaction.plot(x.factor = Nitrogen, 
                                   trace.factor = block, 
                                   response = Yield, 
                                   fun = mean, 
                                   type = "b", 
                                   legend = TRUE,
                                   xlab = "Nitrogen level", 
                                   ylab = "Yield",
                                   pch = c(1, 2, 17, 19)))

## CRD model (ignoring blocks)
potato_CRD.lm <- lm(Yield ~ Nitrogen, data = potato_data)
summary(potato_CRD.lm)
Anova_table(potato_CRD.lm)

## RCBD model (with blocks)
potato_RCBD.lm <- lm(Yield ~ Nitrogen + block, data = potato_data)
summary(potato_RCBD.lm)
Anova_table(potato_RCBD.lm)

## Post-hoc pairwise comparisons (LSD: adjust="none")
emm.CRD  <- emmeans(potato_CRD.lm, specs = "Nitrogen")
pairs(emm.CRD, adjust = "none")

emm.RCBD <- emmeans(potato_RCBD.lm, specs = "Nitrogen")
pairs(emm.RCBD, adjust = "none")

## Plot estimated marginal means
emmip(potato_RCBD.lm, block ~ Nitrogen)
