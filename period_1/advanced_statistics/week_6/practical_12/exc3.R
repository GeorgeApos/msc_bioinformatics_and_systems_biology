# ---------------------------------------------------------
# Exercise 3: Tree volume, height and dbh
# ---------------------------------------------------------

library(car)
library(plotly)
library(tibble)

# Read data
volume_data <- read.table(file = "Tree volume.txt", header = TRUE)
volume_data$clone <- as.factor(volume_data$clone)
str(volume_data)

# Matrix scatterplot
pairs(volume_data[, c(3, 1, 2)])

# Create log(height)
volume_data$lnh <- log(volume_data$height)

# ---------------------------------------------------------
# (a) Fit main regression model
# ---------------------------------------------------------
lnVolume.model1 <- lm(lnvol ~ lndbh + lnh, data = volume_data)
summary(lnVolume.model1)
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
Anova_table(lnVolume.model1)

# Diagnostic plots
plot(lnVolume.model1)

# ---------------------------------------------------------
# (b) Test if β1 = 2 and β2 = 1
# ---------------------------------------------------------
# Example using t-test manually
pt(q = -0.78, df = 82, lower.tail = FALSE)
pt(q = -0.78, df = 82, lower.tail = TRUE)

# Confidence intervals
confint(lnVolume.model1)

# ---------------------------------------------------------
# (c) Confidence interval for given dbh and height
# ---------------------------------------------------------
newdata <- data.frame(lndbh = log(60), lnh = log(25))
predict.lm(object = lnVolume.model1, newdata = newdata,
           se.fit = TRUE, interval = "confidence")

# Back-transform predicted values
exp(8.093)  # fitted value
exp(8.069)  # lower interval
exp(8.117)  # upper interval

# ---------------------------------------------------------
# (d) Include clone effect
# ---------------------------------------------------------
lnVolume.model2 <- lm(lnvol ~ lndbh + lnh + clone, data = volume_data)
summary(lnVolume.model2)
Anova_table(lnVolume.model2)

# ---------------------------------------------------------
# (e) 3D Visualization
# ---------------------------------------------------------
plot_ly(type = "scatter3d", mode = "markers",
        y = volume_data$lnvol,
        x = volume_data$lnh,
        z = volume_data$lndbh,
        color = volume_data$clone)

# ---------------------------------------------------------
# (f) Simpler model: only dbh
# ---------------------------------------------------------
lnVolume.model3 <- lm(lnvol ~ lndbh, data = volume_data)
summary(lnVolume.model3)
