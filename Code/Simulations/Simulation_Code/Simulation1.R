library(tidyverse)
library(tdnnR)

n = 2000

grid <- seq(from = 0, to = 1, by = 0.01)
my_f <- function(x) {
  (5 * x)^2
}
my_F <- Vectorize(FUN = my_f)
f_vals <- my_F(grid)

data_x <- runif(n = n)
data_y <- my_F(data_x) + rnorm(n = n, mean = 0, sd = 2)
data <- cbind(data_y, data_x)

predictions <- unlist(purrr::map(
  .x = grid,
  .f = ~tdnnR::DNN(x = .x, data = data, s = 20)
))

predictions2 <- unlist(purrr::map(
  .x = grid,
  .f = ~tdnnR::TDNN(x = .x, data = data, s1 = 5, s2 = 20)
))
