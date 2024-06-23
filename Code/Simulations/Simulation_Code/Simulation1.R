# This is the first simulation assessing the performance for pointwise estimation

library(tidyverse)
library(tdnnR)
set.seed(1234)

reps <- 10000
n <- 10000
sample_sizes <- c(20, 50, 100, 500, 1000, 5000, 10000)
x <- 0.5

#grid <- seq(from = 0, to = 1, by = 0.01)
my_f <- function(x) {
  (5 * x)^2
}
my_F <- Vectorize(FUN = my_f)
#f_vals <- my_F(grid)

predictions <- matrix(data = NA, nrow = reps, ncol = length(sample_sizes))
predictions2 <- matrix(data = NA, nrow = reps, ncol = length(sample_sizes))

for(i in 1:reps){
  for(j in 1:length(sample_sizes)){
    data_x <- runif(n = sample_sizes[j])
    data_y <- my_F(data_x) + rnorm(n = sample_sizes[j], mean = 0, sd = 2)
    data <- cbind(data_y, data_x)

    predictions[i,j] <- tdnnR::DNN(x = x, data = data[1:sample_sizes[j],], s = 20)
    predictions2[i,j] <- tdnnR::TDNN(x = x, data = data[1:sample_sizes[j],], s1 = 5, s2 = 20)
  }
}

