DNN_Nuisance <- function(x, data, s, asymp_approx_weights = TRUE) {
  # Check inputs for executing function
  point_check(x)
  data_check(data)
  data_point_compatibility_check(x, data[,-2])
  samplingscale_check(s, data)

  # Sort data
  data <- CATE_pre_sort(x = x, data = data)
  untreated_data <- data[data[, 2] == 0, ]
  treated_data <- data[data[, 2] == 1, ]

  W <- as.vector(data[, 2])
  Y0 <- as.vector(data[data[,2] == 0, 1])
  Y1 <- as.vector(data[data[,2] == 1, 1])
  d <- ncol(data) - 2

  n <- length(W)
  n0 <- sum(W == 0)
  n1 <- sum(W == 1)

  # Calculate DNN estimators
  psc_res <- 0
  mu0_res <- 0
  mu1_res <- 0

  psc_factor <- 1
  mu0_factor <- 1
  mu1_factor <- 1

  psc_prefactor <- 0
  mu0_prefactor <- 0
  mu1_prefactor <- 0

  # using exact weights
  if (asymp_approx_weights == FALSE) {
    for (i in 1:(n - s + 1)) {
      index <- n - s + 2 - i
      psc_res <- psc_res + psc_factor * W[index]
      psc_prefactor <- psc_prefactor + psc_factor
      psc_factor <- psc_factor * ((n - index + 1) / i)
    }
    psc_res <- psc_res / psc_prefactor


    for (i in 1:(n0 - s + 1)) {
      index <- n0 - s + 2 - i
      mu0_res <- mu0_res + mu0_factor * Y0[index]
      mu0_prefactor <- mu0_prefactor + mu0_factor
      mu0_factor <- mu0_factor * ((n2 - index + 1) / i)
    }
    mu0_res <- mu0_res / mu0_prefactor

    for (i in 1:(n1 - s + 1)) {
      index <- n1 - s + 2 - i
      mu1_res <- mu1_res + mu1_factor * Y1[index]
      mu1_prefactor <- mu1_prefactor + mu1_factor
      mu1_factor <- mu1_factor * ((n1 - index + 1) / i)
    }
    mu1_res <- mu1_res / mu1_prefactor
  }

  # using approximate weights
  if (asymp_approx_weights == TRUE) {
    psc_alpha <- s / n
    mu0_alpha <- s / n0
    mu1_alpha <- s / n1

    for (i in 1:(n - s + 1)) {
      if (psc_alpha * (1 - psc_alpha)^(i - 1) == 0) {
        break
      }
      psc_res <- psc_res + psc_alpha * (1 - psc_alpha)^(i - 1) * W[i]
      psc_prefactor <- psc_prefactor + psc_alpha * (1 - psc_alpha)^(i - 1)
    }

    psc_res <- psc_res / psc_prefactor

    for (i in 1:(n0 - s + 1)) {
      if (mu0_alpha * (1 - mu0_alpha)^(i - 1) == 0) {
        break
      }
      mu0_res <- mu0_res + mu0_alpha * (1 - mu0_alpha)^(i - 1) * Y0[i]
      mu0_prefactor <- mu0_prefactor + mu0_alpha * (1 - mu0_alpha)^(i - 1)
    }

    mu0_res <- mu0_res / mu0_prefactor

    for (i in 1:(n1 - s + 1)) {
      if (mu1_alpha * (1 - mu1_alpha)^(i - 1) == 0) {
        break
      }
      mu1_res <- mu1_res + mu1_alpha * (1 - mu1_alpha)^(i - 1) * Y1[i]
      mu1_prefactor <- mu1_prefactor + mu1_alpha * (1 - mu1_alpha)^(i - 1)
    }

    mu1_res <- mu1_res / mu1_prefactor
  }

  return(c(psc_res, mu0_res, mu1_res))
}

TDNN_Nuisance <- function(x, data, s1, s2, asymp_approx_weights = TRUE) {
  return(NA)
}
