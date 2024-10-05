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
  Zbase <- gdata$Z_block
  W_block <- gdata$W
  Zim <- gdata$Z.instruments

  K <- length(gdata$cov$Xcols)

  J0 <- nrow(Zbase)
  G <- ncol(Zbase)
  if(!is.null(Zim)) {
    L <- length(Zim$dv)
    G <- G + L # to add residuals
  }
  Td <- length(y_block) # should be nxJ0(J-J0)
  # Initial values
  beta <- matrix(0, K * J0, 1)
  gamma <- rep(0, G)
  sigma2 <- 1
  tau <- rep(1, K)

  # for instruments (if included)
  if(!is.null(gdata$Z.instruments)) {
    sigma.Z <- NULL
    beta.Z <- list()
    for(l in 1:L) {
     sigma.Z <- c(sigma.Z,0)
      beta.Z[[l]] = rep(0,ncol(Zim$Z.im[[l]]))
    }
  }

  # do these outside the main iterations?
  XtW <- Matrix::t(X_block) %*% W_block
  XtWX <- XtW %*% X_block
  XtWy <- XtW %*% y_block

  # Store samples
  beta_samples <- matrix(0, n_iter - burn_in, K * J0)
  gamma_samples <- matrix(0, n_iter - burn_in, G)
  sigma2_samples <- numeric(n_iter - burn_in)
  tau_samples <- matrix(0, n_iter - burn_in, K)

  # for first stage/instruments
  if(!is.null(gdata$Z.instruments)) {
    beta.Z_samples <- sigma.Z_samples <- list()
    for(l in 1:L) {
      sigma.Z_samples[[l]] <- numeric(n_iter - burn_in)
      beta.Z_samples[[l]] <- matrix(0,nrow = n_iter - burn_in,
                                  ncol = ncol(gdata$Z.instruments$Z.im[[l]]))
    }
  }

  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

  for(iter in 1:n_iter) {

    # sample first stage(s) and estimate mean of posterior for data
    if(!is.null(Zim)) {
      gibbsZ <- list()
      .lmat <- NULL
      for(l in 1:L) {
        Y <- gdata$Z_block[,Zim$dv[l]]
        Qz <- Zim$Z.im[[l]]

        # sample
        gibbsZ[[l]] <- gibbs_sampler_one_draw(Y, Qz,
                             beta = beta.Z[[l]],
                             sigma2 = sigma.Z[l])
        # control function: take Z object and replace the 1st stage dv
        # note: we only need to do this once for each endogenous variable
        # regardless of whether it is subsequently squared or interacted with
        # another variable. However, it must include correspondingly more
        # instruments

        # samples
        sigma.Z[l] <- gibbsZ[[l]]$sigmaZ2
        beta.Z[[l]] <- gibbsZ[[l]]$betaZ

        # residuals appended to covariates in 2nd stage of hierarchical model
        .l <- gibbsZ[[l]]$residuals
        colnames(.l) <- paste0("lres",l)
        .lmat <- cbind(.lmat, .l)

      }
      Z <- cbind(Zbase, .lmat)
    } else Z <- Zbase

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

      if(!is.null(gdata$Z.instruments)) {
        for(l in 1:L) {
        beta.Z_samples[[l]][iter-burn_in,] <- beta.Z[[l]]
        sigma.Z_samples[[l]][iter-burn_in] <- sigma.Z[l]
        }
      }

    }
    # Update progress bar
    setTxtProgressBar(pb, iter)
  }

  # this helps for later interpretation
  colnames(gamma_samples) <- colnames(Z)

  if(!is.null(gdata$Z.instruments)) {

    output <- list(beta_samples = beta_samples, gamma_samples = gamma_samples,
                   sigma2_samples = sigma2_samples, tau_samples = tau_samples,
                   beta.Z_samples = beta.Z_samples,
                   sigma.Z_samples = sigma.Z_samples)
  } else {
    output <- list(beta_samples = beta_samples, gamma_samples = gamma_samples,
                   sigma2_samples = sigma2_samples, tau_samples = tau_samples)
  }
  return(output)
}


#' Perform one iteration of the Gibbs sampler for an OLS model
#'
#' This function performs one iteration of a Gibbs sampler for a simple OLS
#' (Ordinary Least Squares) regression model. It samples the coefficients
#' (beta) and the variance of the errors (sigma^2) given the current data and
#' parameters.
#'
#' @param y A vector of the response variable values (n x 1).
#' @param X A matrix of the explanatory variables (n x p), including the intercept.
#' @param beta A vector of current values of the regression coefficients (p x 1).
#' @param sigma2 The current value of the variance of the errors.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{beta_sample}{A vector of the sampled regression coefficients.}
#'   \item{sigma2_sample}{The sampled variance of the errors.}
#'   \item{residuals}{The residuals from the regression, i.e., y - X * beta.}
#' }
#'
#' @examples
#' # Simulated data example
#' set.seed(123)
#' n <- 100
#' p <- 2
#' X <- cbind(1, rnorm(n))  # Design matrix (including intercept term)
#' beta_true <- c(2, 3)  # True beta coefficients
#' sigma2_true <- 1  # True variance of the error
#'
#' # Generate y according to the model y = X * beta + noise
#' y <- X %*% beta_true + rnorm(n, mean = 0, sd = sqrt(sigma2_true))
#'
#' # Initial values for beta and sigma^2
#' beta_init <- rep(0, p)
#' sigma2_init <- 1
#'
#' # Run one iteration of the Gibbs sampler
#' sample_result <- gibbs_sampler_one_draw(y, X, beta_init, sigma2_init)
#'
#' # Output
#' print(sample_result$beta_sample)
#' print(sample_result$sigma2_sample)
#' print(sample_result$residuals[1:10])
#' @importFrom MASS mvrnorm
#'
#' @export
gibbs_sampler_one_draw <- function(y, X, beta, sigma2) {
  n <- nrow(X)
  p <- ncol(X)

  # Precompute some matrix operations for efficiency
  XtX_inv <- solve(t(X) %*% X)

  # 1. Sample beta | y, X, sigma^2
  beta_mean <- XtX_inv %*% t(X) %*% y
  beta_var <- sigma2 * XtX_inv
  beta_sample <- mvrnorm(1, beta_mean, beta_var)  # Draw beta from the multivariate normal distribution

  # 2. Sample sigma^2 | y, X, beta
  residuals <- y - X %*% beta_sample
  alpha <- (n / 2)
  beta_param <- sum(residuals^2) / 2
  sigma2_sample <- 1 / rgamma(1, shape = alpha, rate = beta_param)  # Draw sigma^2 from the inverse-gamma

  # Return beta sample, sigma^2 sample, and residuals
  return(list(betaZ = beta_sample,
              sigmaZ2 = sigma2_sample,
              residuals = residuals))
}

