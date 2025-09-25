## Read data into an object named firstGrass_data
# complete the code
# Note that sep="\t" indicates that columns are separated by tabs
# reading in a header can be TRUE or FALSE
firstGrass_data <- read.table(file = "First cut of Grass.txt", header = TRUE )

# Create the most appropriate plot (1, 2 or 3)
# option 1 boxplot
with(firstGrass_data, boxplot(Yield ~ Mix))

# option 2 jitterplot

ggplot(firstGrass_data, aes(x=Mix, y=Yield)) + 
  geom_jitter(position=position_jitter(0.1))

#label: create dummies

# Create dummies for the levels of 'Mix"
firstGrass_data <- dummy_cols(firstGrass_data, select_columns ="Mix" )
head(firstGrass_data, 15)

