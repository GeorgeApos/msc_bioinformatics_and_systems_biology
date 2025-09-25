library(doBy)
library(emmeans)
library(ggplot2)

# run Anova_table type 2 function
# Advanced statistics CP 5 thru 12
# This function extends the anova table that is produced
# by Anova() from the car library
# It is set to produce type II sums of squares.
# argument model: an object as produced by the lm() function

library(car)
library(tibble)
Anova_table <- function(model) {
  
  #get the TypeII Anova table 
  At2<-Anova(model, Type =2)
  
  #construct SStot and dftot from a TypeI Anova table
  At1<-anova(model)
  SStotal <- sum(At1[,2])
  dftotal <- sum(At1[,1])
  Total<-as.data.frame(cbind(SStotal, dftotal,"",""))
  rownames(Total)<-"Total"
  names(Total) <- names(At2)
  
  Avtable <- rbind(At2,Total)  
  
  for (c in 1:4) {
    Avtable[,c]<-as.numeric(Avtable[,c])
  }  
  rows<-dim(Avtable)[1] 
  MS<-vector(length=rows)
  MS<- as.numeric(Avtable[,1])/as.numeric(Avtable[,2])
  
  Avtable <- add_column(Avtable,MS,.after=2)
  
  return(Avtable) 
}

# read in the data 
hostility_data <- read.table(file="Hostility.txt", header = TRUE)
head(hostility_data)

hostility_data$method <- as.factor(hostility_data$method) # if necessary

# get summary statistics for each method separately using doBy
summaryBy(Hscore ~ method, data = hostility_data,
          FUN = function(x) { c(m = mean(x), s = sd(x), se = sd(x)/sqrt(length(x)), n=length(x) )  } )
## estimate the ANOVA-model with lm(); complete the code
hostility.lm <- lm(formula = Hscore ~ method, data = hostility_data  )

## get output of the models using summary
summary(hostility.lm)

## get output of the models using anova
Anova_table(hostility.lm)

## use emmeans library to estimate expected marginal means
# Install once
install.packages("emmeans")

# Load in each session
library(emmeans)

# request EM-means for the response variable per method
emmeans(hostility.lm, specs = "method")

# plot the expected marginal means with interval
plot(emmeans(hostility.lm, specs = "method"), horizontal = FALSE)

# compute pairwise differences
# protected LSD : adjust = "none", so no multiple comparison adjustment
# default = "tukey"
pairs(emmeans(hostility.lm, specs = "method"), adjust = "none")

pairs(emmeans(hostility.lm, specs = "method"), adjust = "tukey")

