#' Gibbs Sampling Function
#'
#' This function runs a Gibbs sampler for a given set of parameters.
#'
#' @param gdata An object produced by augsynth "multisynth"
#'    gdata$y A numeric vector of responses.
#'    gdata$X A design matrix.
#'    gdata$W A weight matrix.
#'    gdata$Z A covariate matrix.
#' @param n_iter Number of iterations to run.
#' @param burn_in Number of iterations to discard as burn-in.
#' @return A list of MCMC samples.
#' @importFrom stats rgamma rnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
gibbs_sampling <- function(gdata,
                           n_iter = 1000,
                           burn_in = 500) {

  y_block <- gdata$Y_block
  X_block <- gdata$X_block
  Z_block <- gdata$Z_block
  W_block <- gdata$W
  Z <- gdata$Z

  K <- length(gdata$cov$Xcols)

  J0 <- nrow(Z)
  G <- ncol(Z)
  Td <- length(y_block) # should be nxJ0(J-J0)

  # Initial values
  beta <- matrix(0, K * J0, 1)
  gamma <- rep(0, G)
  sigma2 <- 1
  tau <- rep(1, K)

  # do these outside the main iterations?
  XtW <- Matrix::t(X_block) %*% W_block
  XtWX <- XtW %*% X_block
  XtWy <- XtW %*% y_block

  # Store samples
  beta_samples <- matrix(0, n_iter - burn_in, K * J0)
  gamma_samples <- matrix(0, n_iter - burn_in, G)
  sigma2_samples <- numeric(n_iter - burn_in)
  tau_samples <- matrix(0, n_iter - burn_in, K)

  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

  for(iter in 1:n_iter) {
    # set up Zstar
    dZ <- rep(0,K)
    dZ[gdata$cov$intX] <- 1
    Z_star <- Matrix::kronecker(Z, dZ)

    # sample beta given y, X, W, sigma2, gamma, tau
    Sigma_beta_prior_inv <- Matrix::kronecker(diag(J0), diag(1/tau)) # prior inverse of beta_sigma
    q <- Z_star %*% gamma
    A <- XtWX / sigma2 + Sigma_beta_prior_inv
    b <- XtWy / sigma2 + Sigma_beta_prior_inv %*% q
    cA <- Matrix::chol(A)
    A_inv <- Matrix::chol2inv(cA)
    mu_beta <- A_inv%*%b

    beta <- rMVNormCovariance(1, mu = mu_beta, Sigma = A_inv)

    beta_matrix <- matrix(beta, ncol = K, byrow = TRUE)

    # Sample gamma given beta and tau - using ivGibbs
    beta_2 <- beta_matrix[, gdata$cov$intX]
    V_gamma <- Matrix::solve(t(Z) %*% Z / tau[2]^2 + diag(G))
    m_gamma <- V_gamma %*% (t(Z) %*% beta_2 / tau[2]^2)
    gamma <- rMVNormCovariance(1, mu = m_gamma, Sigma = V_gamma)

    # Sample sigma2 given y and beta
    residuals <- y_block - X_block %*% beta
    alpha <- length(y_block) / 2
    beta_param <- sum((Matrix::t(residuals) %*% W_block %*% residuals)) / 2
    sigma2 <- 1 / rgamma(1, shape = alpha, rate = beta_param)

    # Sample tau given beta and gamma
    for (k in 1:K) {
      tau[k] <- sqrt(1 / rgamma(1, shape = (J0 / 2),
                                rate = (sum((beta_matrix[, k] - Z_star[((1:J0) - 1) * K + k, ] %*% gamma)^2) / 2)))
    }

    # Store samples after burn-in

    if(iter == burn_in) cat("\nBeginning sampling after burnin.")
    if (iter > burn_in) {
      beta_samples[iter - burn_in, ] <- as.numeric(beta)
      gamma_samples[iter - burn_in, ] <- as.numeric(gamma)
      sigma2_samples[iter - burn_in] <- as.numeric(sigma2)
      tau_samples[iter - burn_in, ] <- as.numeric(tau)
    }
    # Update progress bar
    setTxtProgressBar(pb, iter)
  }

  return(list(beta_samples = beta_samples, gamma_samples = gamma_samples,
              sigma2_samples = sigma2_samples, tau_samples = tau_samples))
}
