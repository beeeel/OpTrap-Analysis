# First attempt at data analysis in R
# load the data

# Need this for read_tsv
library(readr) 
# It's easier to save in MATLAB as tab separated.
data <- read_tsv("~/Data/phd/OpTrap/accumulator_2021_07_21-latB/sim_a07D0025.csv", col_names=FALSE)

# Get the data out of the tibble into vectors
data$X1 -> tau
data$X2 -> msd
