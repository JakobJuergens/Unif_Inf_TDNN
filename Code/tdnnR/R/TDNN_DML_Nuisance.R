#' This function estimates the propensity score, mu0 and mu1
#' at a given point
#'
#' @param x Point of interest
#' @param data Data matrix for the CATE setting
#' @param s1 Subsampling Scale 1
#' @param s2 Subsampling Scale 2
#' @param asymp_approx_weights True or False whether to use an asymptotic
#' approximation for the weights (default value = FALSE)
#'
TDNN_Nuisance <- function(x, data, s1, s2, asymp_approx_weights = TRUE) {

  # Check that s2 > s1 and reverse if not
  if (s1 > s2) {
    tmp <- s1
    s1 <- s2
    s2 <- tmp
  }

  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data[,-2])
  samplingscale_check(s, data)

  # Sort data
  data <- CATE_pre_sort(x = x, data = data)$data
  untreated_data <- data[data[, 2] == 0, ]
  treated_data <- data[data[, 2] == 1, ]

  W <- as.vector(data[, 2])
  Y0 <- as.vector(data[data[,2] == 0, 1])
  Y1 <- as.vector(data[data[,2] == 1, 1])
  d <- ncol(data) - 2

  n <- length(W)
  n0 <- sum(W == 0)
  n1 <- sum(W == 1)

  # Calculate weights for two-scale DNN
  w1 <- (1 - (s1 / s2)^(-2 / d))^(-1)
  w2 <- 1 - w1

  # Calculate TDNN estimators
  psc_res1 <- 0
  mu0_res1 <- 0
  mu1_res1 <- 0
  psc_res2 <- 0
  mu0_res2 <- 0
  mu1_res2 <- 0

  psc_factor1 <- 1
  mu0_factor1 <- 1
  mu1_factor1 <- 1
  psc_factor2 <- 1
  mu0_factor2 <- 1
  mu1_factor2 <- 1

  psc_prefactor1 <- 0
  mu0_prefactor1 <- 0
  mu1_prefactor1 <- 0
  psc_prefactor2 <- 0
  mu0_prefactor2 <- 0
  mu1_prefactor2 <- 0

  # using exact weights
  if(asymp_approx_weights == FALSE){

    for (i in 1:(s2 - s1)) {
      psc_index <- n - s1 + 2 - i
      mu0_index <- n0 - s1 + 2 - i
      mu1_index <- n1 - s1 + 2 - i

      psc_res1 <- psc_res1 + psc_factor1 * W[psc_index]
      psc_prefactor1 <- psc_prefactor1 + psc_factor1
      psc_factor1 <- psc_factor1 * ((n - psc_index + 1) / i)

      mu0_res1 <- mu0_res1 + mu0_factor1 * Y0[mu0_index]
      mu0_prefactor1 <- mu0_prefactor1 + mu0_factor1
      mu0_factor1 <- mu0_factor1 * ((n - mu0_index + 1) / i)

      mu1_res1 <- mu1_res1 + mu1_factor1 * Y1[mu1_index]
      mu1_prefactor1 <- mu1_prefactor1 + mu1_factor1
      mu1_factor1 <- mu1_factor1 * ((n - mu1_index + 1) / i)
    }

    for (i in 1:(n - s2 + 1)) {
      psc_index <- n - s2 + 2 - i

      psc_res1 <- psc_res1 + psc_factor1 * W[psc_index]
      psc_res2 <- psc_res2 + psc_factor2 * W[psc_index]
      psc_prefactor1 <- psc_prefactor1 + psc_factor1
      psc_prefactor2 <- psc_prefactor2 + psc_factor2

      psc_factor1 <- psc_factor1 * ((n - psc_index + 1) / (i + s2 - s1))
      psc_factor2 <- psc_factor2 * ((n - psc_index + 1) / i)
    }

    for (i in 1:(n0 - s2 + 1)) {
      mu0_index <- n - s2 + 2 - i

      mu0_res1 <- mu0_res1 + mu0_factor1 * Y0[mu0_index]
      mu0_res2 <- mu0_res2 + mu0_factor2 * Y0[mu0_index]
      mu0_prefactor1 <- mu0_prefactor1 + mu0_factor1
      mu0_prefactor2 <- mu0_prefactor2 + mu0_factor2

      mu0_factor1 <- mu0_factor1 * ((n - mu0_index + 1) / (i + s2 - s1))
      mu0_factor2 <- mu0_factor2 * ((n - mu0_index + 1) / i)
    }

    for (i in 1:(n1 - s2 + 1)) {
      mu1_index <- n - s2 + 2 - i

      mu1_res1 <- mu1_res1 + mu1_factor1 * Y1[mu1_index]
      mu1_res2 <- mu1_res2 + mu1_factor2 * Y1[mu1_index]
      mu1_prefactor1 <- mu1_prefactor1 + mu1_factor1
      mu1_prefactor2 <- mu1_prefactor2 + mu1_factor2

      mu1_factor1 <- mu1_factor1 * ((n - mu1_index + 1) / (i + s2 - s1))
      mu1_factor2 <- mu1_factor2 * ((n - mu1_index + 1) / i)
    }
  }

  # using asymptotically approximated weights
  if(asymp_approx_weights == TRUE){
    psc_alpha_1 <- s1/n
    psc_alpha_2 <- s2/n
    mu0_alpha_1 <- s1/n0
    mu0_alpha_2 <- s2/n0
    mu1_alpha_1 <- s1/n1
    mu1_alpha_2 <- s2/n1

    # Propensity Score
    for (i in 1:(n - s1 + 1)) {
      if(psc_alpha_1*(1-psc_alpha_1)^(i-1) == 0){break}
      psc_res1 <- psc_res1 + psc_alpha_1*(1-psc_alpha_1)^(i-1) * W[i]
      psc_prefactor1 <- psc_prefactor1 + psc_alpha_1*(1-psc_alpha_1)^(i-1)
    }

    for (i in 1:(n - s2 + 1)) {
      if(psc_alpha_2*(1-psc_alpha_2)^(i-1) == 0){break}
      psc_res2 <- psc_res2 + psc_alpha_2*(1-psc_alpha_2)^(i-1) * W[i]
      psc_prefactor2 <- psc_prefactor2 + psc_alpha_2*(1-psc_alpha_2)^(i-1)
    }

    # mu_0
    for (i in 1:(n0 - s1 + 1)) {
      if(mu0_alpha_1*(1-mu0_alpha_1)^(i-1) == 0){break}
      mu0_res1 <- mu0_res1 + mu0_alpha_1*(1-mu0_alpha_1)^(i-1) * Y0[i]
      mu0_prefactor1 <- mu0_prefactor1 + mu0_alpha_1*(1-mu0_alpha_1)^(i-1)
    }

    for (i in 1:(n0 - s2 + 1)) {
      if(mu0_alpha_2*(1-mu0_alpha_2)^(i-1) == 0){break}
      mu0_res2 <- mu0_res2 + mu0_alpha_2*(1-mu0_alpha_2)^(i-1) * Y0[i]
      mu0_prefactor2 <- mu0_prefactor2 + mu0_alpha_2*(1-mu0_alpha_2)^(i-1)
    }

    # mu_1
    for (i in 1:(n1 - s1 + 1)) {
      if(mu1_alpha_1*(1-mu1_alpha_1)^(i-1) == 0){break}
      mu1_res1 <- mu1_res1 + mu1_alpha_1*(1-mu1_alpha_1)^(i-1) * Y1[i]
      mu1_prefactor1 <- mu1_prefactor1 + mu1_alpha_1*(1-mu1_alpha_1)^(i-1)
    }

    for (i in 1:(n1 - s2 + 1)) {
      if(mu1_alpha_2*(1-mu1_alpha_2)^(i-1) == 0){break}
      mu1_res2 <- mu1_res2 + mu1_alpha_2*(1-mu1_alpha_2)^(i-1) * Y1[i]
      mu1_prefactor2 <- mu1_prefactor2 + mu1_alpha_2*(1-mu1_alpha_2)^(i-1)
    }
  }

  psc_res1 <- psc_res1/psc_prefactor1
  psc_res2 <- psc_res2/psc_prefactor2
  mu0_res1 <- mu0_res1/mu0_prefactor1
  mu0_res2 <- mu0_res2/mu0_prefactor2
  mu1_res1 <- mu1_res1/mu1_prefactor1
  mu1_res2 <- mu1_res2/mu1_prefactor2

  # Combine results
  psc_res <- w1 * psc_res1 + w2 * psc_res2
  mu0_res <- w1 * mu0_res1 + w2 * mu0_res2
  mu1_res <- w1 * mu1_res1 + w2 * mu1_res2

  # Careful! This currently returns a DNN estimate for the Propensity Score
  return(c(psc_res1, mu0_res, mu1_res))
}
