# This is the first simulation assessing the performance for pointwise estimation
# It replicates setting 1 from Demirkaya et al. (2024)

library(tidyverse)
library(tdnnR)
library(MASS)

set.seed(1234)

reps <- 10000
sample_sizes <- c(50, 100, 500, 1000, 5000, 10000)
test_point <- c(0.5, -0.5, 0.5)

my_f <- function(i, data) {
  return((data[i, 1]-1)^2 + (data[i,2]+1)^3 - 3*data[i,3])
}
my_F <- Vectorize(FUN = my_f, vectorize.args = list('i'))

predictions <- matrix(data = NA, nrow = reps, ncol = length(sample_sizes))
predictions2 <- matrix(data = NA, nrow = reps, ncol = length(sample_sizes))

for (i in 1:reps) {
  if (i %% 10 == 0) {
    print(paste0("Repetition: ", i))
  }
  for (j in 1:length(sample_sizes)) {
    data_x <- MASS::mvrnorm(n = sample_sizes[j], mu = c(0,0,0), Sigma = diag(c(1,1,1)))
    data_y <- my_F(data = data_x, i = 1:sample_sizes[j]) + rnorm(n = sample_sizes[j], mean = 0, sd = 1)
    data <- cbind(data_y, data_x)

    predictions[i, j] <- tdnnR::DNN(x = test_point, data = data, s = 20, standardize = FALSE)
    predictions2[i, j] <- tdnnR::TDNN(x = test_point, data = data, s1 = 20, s2 = 50, standardize = FALSE)
  }
}

pred_tibble <- as_tibble(predictions)
names(pred_tibble) <- unlist(purrr::map(
  .x = 1:length(sample_sizes),
  .f = ~ paste0("n = ", sample_sizes[.x])
))
pred_tibble <- pred_tibble |> pivot_longer(cols = everything(), names_to = 'Sample Size')
pred_tibble$Type <- 'DNN'
pred_tibble2 <- as_tibble(predictions2)
names(pred_tibble2) <- unlist(purrr::map(
  .x = 1:length(sample_sizes),
  .f = ~ paste0("n = ", sample_sizes[.x])
))
pred_tibble2 <- pred_tibble2 |> pivot_longer(cols = everything(), names_to = 'Sample Size')
pred_tibble2$Type <- 'TDNN'
pred_tibble_comb <- rbind(pred_tibble, pred_tibble2)
pred_tibble_comb$`Sample Size` <- factor(pred_tibble_comb$`Sample Size`,
                                         levels = unique(pred_tibble_comb$`Sample Size`))

TDNN_DNN_plot <- ggplot(data = pred_tibble_comb) +
  geom_violin(aes(x = `Sample Size`, y = value, fill = `Type`),
              trim = FALSE) +
  ylab("Estimate") +
  #labs(title = '(T)DNN Estimator for s = 20, s1 = 20, s2 = 50',
  #     subtitle = 'for different sample sizes') +
  geom_hline(yintercept = my_f(data = t(as.matrix(test_point)), i = 1), linetype = 2) +
  scale_fill_brewer(palette="Dark2") +
  theme_minimal() +
  theme(legend.position="bottom")

