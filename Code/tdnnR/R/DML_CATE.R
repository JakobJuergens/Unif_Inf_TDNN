DNN_DML2 <- function(x, data, s, n_folds,
                     asymp_approx_weight = TRUE) {
  # Create List of Fold indices to Estimate Propensity Scores and Regression Functions
  folds <- generate_folds(n_obs = nrow(data), n_folds = n_folds)
  nuisance_par_ests <- matrix(data = NA, nrow = nrow(data), ncol = 3)

  # Identify treated and untreated units
  untreated <- which(data[,2] == 0)
  treated <- which(data[,2] == 1)

  # Estimate Nuisance Parameters
  for (i in 1:n_folds) {
    tmp <- purrr::map(
        .x = folds[[i]],
        .f = ~ DNN_Nuisance(x = data[.x, -(1:2)], data = data[-folds[[i]], ], s = s,
                            asymp_approx_weights = asymp_approx_weights)
      )
    for(j in 1:length(folds[[i]])){
      nuisance_par_ests[folds[[i]][j], ] <- tmp[[j]]
    }
  }

  # Estimate Parameter of Interest
}

TDNN_DML2 <- function(x, data, s1, s2, n_folds) {
  n_obs <- nrow(data)
  # Create Folds to Estimate Propensity Scores
}
