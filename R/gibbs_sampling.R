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
#' @param Z_cov_dense logical - whether use dense covariance
#' @param run_selection_gibbs logical - whether to run selection eqn
#' @return A list of MCMC samples.
#' @importFrom stats rgamma rnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix kronecker as.matrix
#' @export
gibbs_sampling <- function(gdata,
                           n_iter = 1000,
                           burn_in = 500,
                           Z_cov_dense = FALSE,
                           run_selection_gibbs = FALSE) {

  y_block <- gdata$Y_block
  X_block <- gdata$X_block
  Zbase <- gdata$Z_block
  W_block <- gdata$W
  Zim <- gdata$Z.instruments

  # design matrix for probit selection equation
  GRX <- gdata$GRX

  K <- length(gdata$cov$Xcols)

  J0 <- gdata$cov$J0
  J <- gdata$cov$J
  G <- ncol(Zbase)
  if(!is.null(Zim)) {
    L <- length(Zim$dv)
    G <- G + L # to add residuals
  }
  Td <- length(y_block) # should be nxJ0(J-J0)

  #### Set up Initial values
  beta <- matrix(0, K * J0, 1)
  gamma <- rep(0, G)
  sigma2 <- 1
  tau <- rep(1, K)

  # for selection equation
  if(run_selection_gibbs) {
    beta_fs <- matrix(0, ncol(GRX), 1)
    sigma2_fs <- 1
  }

  # for instruments for second stage (if included)
  if(!is.null(Zim)) {
    sigma.Z <- list()
    if(Z_cov_dense) {
      sigma.Z <- initialize_dense_covariance(L)
    } else {
      for(l in 1:L) {
        # one for each
        sigma.Z <- c(sigma.Z,1)
      }
      sigma.Z <- diag(sigma.Z)
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
  if(run_selection_gibbs) {
    beta_fs_samples <- matrix(0, nrow = n_iter - burn_in, ncol = ncol(GRX))
    sigma2_fs_samples <- numeric(n_iter - burn_in)
  }

  # for first stage/instruments
  if(!is.null(Zim)) {
    beta.Z_samples <- list()
    if(!Z_cov_dense) {
      sigma.Z_samples <- matrix(0, nrow= n_iter - burn_in, ncol = L)
    } else sigma.Z_samples <- array(0, dim = c(L,L,n_iter - burn_in))

    for(l in 1:L) beta.Z_samples[[l]] <- matrix(0,nrow = n_iter - burn_in,
                                  ncol = ncol(Zim$Z.im[[l]]))
    }


  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

  for(iter in 1:n_iter) {

    ############################ Start selection block
    if(run_selection_gibbs) {


      # first stage selection equation (this needs a Gibbs step function)
      .bd <- gibbs_binomial_probit(GRX, gdata$GRY,
                                      beta = beta_fs,
                                      Sigma_prior = diag(ncol(GRX)),
                                      mu_prior = rep(0,ncol(GRX)),
                                      sigma2 = sigma2_fs,
                                      alpha_0 = 2,
                                      beta_0 = 2)

      beta_fs <- .bd$beta
      sigma2_fs <- .bd$sigma2

      resid <- gdata$GRY - pnorm( GRX %*% beta_fs)
      residfake <- rnorm(nrow(X_block))

      # this step is slow and requires a lot of memory
      X_block <- replace_kth_column_block_diagonal_fast(X_block,
                                                        residfake,
                                                        T = gdata$cov$T * (1+J-J0),
                                                        K = K,
                                                        N = J0,
                                                        k = which(gdata$cov$Xcols=="GR"))
      # and this part a bit slow
      XtW <- Matrix::t(X_block) %*% W_block
      XtWX <- XtW %*% X_block
      XtWy <- XtW %*% y_block

    }
    ########################### End selection block

    ########################### Start endogenous regression first stage
    # sample first stage(s) and estimate mean of posterior for data
      if(!is.null(Zim)) {
        gibbsZ <- Ylist <- Qzlist <- list()
      .lmat <- NULL
      # assemble the multivariate endogenous variables list
      for(l in 1:L) {
        Ylist[[l]] <- gdata$Z_block[,Zim$dv[l]]
        Qzlist[[l]] <- Zim$Z.im[[l]]
      }
        # sample

        if(!Z_cov_dense) {
          gibbsZ <- gibbs_sampler_one_draw_cov(Ylist, Qzlist,
                                               sigma_param = sigma.Z,
                                               covariance_type = "diagonal")
          } else gibbsZ <- gibbs_sampler_one_draw_cov(Ylist, Qzlist,
                             sigma_param = sigma.Z,
                             covariance_type = "dense")

        # control function: take Z object and replace the 1st stage dv
        # note: we only need to do this once for each endogenous variable
        # regardless of whether it is subsequently squared or interacted with
        # another variable. However, it must include correspondingly more
        # instruments

        # samples - cov is either list if diagonal or matrix if dense
        sigma.Z <- gibbsZ$sigma2_samples
        beta.Z <- gibbsZ$beta_samples

        # residuals appended to covariates in 2nd stage of hierarchical model
        .lmat <- gibbsZ$residuals
        colnames(.lmat) <- paste0("lres",1:L)

        Z <- Matrix::as.matrix(cbind(Zbase, .lmat))
      } else Z <- Zbase

    ######################### end endogenous regressions first stage

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

      if(!is.null(Zim)) {
        for(l in 1:L) beta.Z_samples[[l]][iter-burn_in,] <- beta.Z[[l]]
        if(Z_cov_dense) {
          sigma.Z_samples[,,iter-burn_in] <- sigma.Z
        } else sigma.Z_samples[iter-burn_in,] <- diag(sigma.Z)
      }

      if(run_selection_gibbs){
        beta_fs_samples[iter-burn_in,] <- beta_fs
        sigma2_fs_samples[iter-burn_in] <- sigma2_fs
      }
    }

    # Update progress bar
    setTxtProgressBar(pb, iter)
  }

  # this helps for later interpretation
  colnames(gamma_samples) <- colnames(Z)

  if(!is.null(Zim)) {
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
#' parameters. Note there is currently no prior that can be passed here.
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


#' Gibbs Sampler for Multivariate Regression with Covariance Option
#'
#' This function performs a single Gibbs sampling draw for a multivariate regression model with an option
#' to use either a dense or diagonal covariance matrix for the outcomes. Each outcome has its own set of predictors.
#'
#' The model assumes that each outcome \eqn{y_j} follows a regression model:
#' \deqn{y_j = X_j \beta_j + \epsilon_j, \quad \epsilon_j \sim \mathcal{N}(0, \Sigma)}
#' where \eqn{X_j} is the predictor matrix for outcome \eqn{y_j}, \eqn{\beta_j} is the vector of regression
#' coefficients, and \eqn{\Sigma} is the covariance matrix of the residuals. The covariance matrix can either
#' be diagonal (independent outcomes) or dense (allowing correlation between outcomes).
#'
#' @param y_list A list of outcome vectors. Each element \eqn{y_j} in the list corresponds to the outcomes for one dependent variable.
#' @param X_list A list of predictor matrices. Each element \eqn{X_j} in the list corresponds to the predictor matrix for the corresponding outcome \eqn{y_j}.
#' @param sigma_param A list of initial values for the residual variances \eqn{\sigma^2_j} for each outcome when `covariance_type = "diagonal"`, or a dense covariance matrix when `covariance_type = "dense"`.
#' @param covariance_type A string specifying the type of covariance matrix to use. Can be either `"diagonal"` (no correlation between outcomes) or `"dense"` (allows correlation between outcomes).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{beta_samples}}{A list of sampled regression coefficients \eqn{\beta_j} for each outcome.}
#'   \item{\code{sigma2_samples}}{A list of sampled residual variances \eqn{\sigma^2_j} for each outcome (diagonal case), or the sampled covariance matrix \eqn{\Sigma} (dense case).}
#'   \item{\code{residuals}}{A list of residuals for each outcome \eqn{r_j = y_j - X_j \beta_j}.}
#' }
#'
#' @examples
#' # Simulate data for two outcomes with different predictor matrices
#' set.seed(123)
#' n <- 100  # Number of observations
#' p <- 2    # Number of outcomes
#'
#' # Create predictor matrices X for each outcome
#' X1 <- cbind(1, rnorm(n))
#' X2 <- cbind(1, rnorm(n), rnorm(n))
#'
#' # True beta values
#' beta1 <- c(1, 0.5)
#' beta2 <- c(1, -0.3, 0.2)
#'
#' # Simulate outcomes
#' y1 <- X1 %*% beta1 + rnorm(n)
#' y2 <- X2 %*% beta2 + rnorm(n)
#'
#' # Prepare data for the Gibbs sampler
#' y_list <- list(y1, y2)
#' X_list <- list(X1, X2)
#' sigma_param <- list(1, 1)  # Diagonal case: Initial values for sigma^2
#'
#' # Run the Gibbs sampler with diagonal covariance
#' result_diag <- gibbs_sampler_one_draw_cov(y_list, X_list, sigma_param, covariance_type = "diagonal")
#' print(result_diag)
#'
#' # Simulate a dense covariance matrix
#' sigma_param_dense <- matrix(c(1, 0.5, 0.5, 1), 2, 2)  # Dense covariance matrix
#'
#' # Run the Gibbs sampler with dense covariance
#' result_dense <- gibbs_sampler_one_draw_cov(y_list, X_list, sigma_param_dense, covariance_type = "dense")
#' print(result_dense)
#'
#' @export
gibbs_sampler_one_draw_cov <- function(y_list, X_list, sigma_param,
                                       covariance_type = "diagonal") {
  p <- length(y_list)  # Number of outcomes
  n <- nrow(X_list[[1]])  # Number of observations (assumed to be the same for all outcomes)

  # Storage for results
  beta_samples <- vector("list", p)
  residuals_matrix <- matrix(0, nrow = n, ncol = p)  # Store residuals for all outcomes

  if (covariance_type == "diagonal") {
    # Diagonal case: sigma_param is a list of sigma^2 values
    sigma2_list <- sigma_param
  } else if (covariance_type == "dense") {
    # Dense case: sigma_param is the dense covariance matrix
    Sigma <- sigma_param
    Sigma_inv <- solve(Sigma)  # Inverse of the covariance matrix
  } else {
    stop("Invalid covariance type. Use 'diagonal' or 'dense'.")
  }

  # Step 1: Sample beta and sigma^2 for each outcome (or use dense covariance)
  for (j in 1:p) {
    X_j <- X_list[[j]]
    y_j <- y_list[[j]]

    # Precompute some matrix operations for efficiency
    XtX_inv_j <- solve(t(X_j) %*% X_j)

    # 1. Sample beta_j | y_j, X_j, sigma^2_j or Sigma (if dense)
    if (covariance_type == "diagonal") {
      beta_mean_j <- XtX_inv_j %*% t(X_j) %*% y_j
      beta_var_j <- sigma2_list[[j]] * XtX_inv_j
      beta_samples[[j]] <- MASS::mvrnorm(1, beta_mean_j, beta_var_j)  # Draw beta from the multivariate normal distribution
    } else if (covariance_type == "dense") {
      beta_mean_j <- XtX_inv_j %*% t(X_j) %*% y_j
      beta_var_j <- Sigma_inv[j, j] * XtX_inv_j  # Use the inverse of Sigma for variance
      beta_samples[[j]] <- MASS::mvrnorm(1, beta_mean_j, beta_var_j)  # Draw beta from the multivariate normal distribution
    }

    # 2. Compute residuals and store them
    residuals_j <- y_j - X_j %*% beta_samples[[j]]
    residuals_matrix[, j] <- residuals_j  # Store residuals for this outcome

    # 3. Sample sigma^2_j | y_j, X_j, beta_j (diagonal case)
    if (covariance_type == "diagonal") {
      alpha <- n / 2
      beta_param <- sum(residuals_j^2) / 2
      sigma2_list[j] <- 1 / rgamma(1, shape = alpha, rate = beta_param)  # Sample sigma^2 from inverse-gamma
    }
  }

  # Step 4: If dense covariance, sample the covariance matrix
  if (covariance_type == "dense") {
    # Degrees of freedom for inverse-Wishart
    nu <- n + p + 1
    # Scale matrix for inverse-Wishart
    S <- t(residuals_matrix) %*% residuals_matrix
    sigma_sample <- MCMCpack::riwish(nu, S)  # Sample covariance matrix from inverse-Wishart
  } else {
    sigma_sample <- sigma2_list  # For diagonal case, return the sampled sigma^2 values from diagonal
  }

  # Return results: beta samples, residuals, and covariance matrix (if dense)
  return(list(beta_samples = beta_samples,
              sigma2_samples = sigma_sample,
              residuals = as.data.frame(residuals_matrix),
              covariance = if (covariance_type == "dense") sigma_sample else NULL))
}



#' Gibbs Sampling for a Binomial Probit Model with sigma^2
#'
#' This function performs one Gibbs sampling iteration for a binomial probit model,
#' with an additional step to update the variance parameter \eqn{\sigma^2}.
#' The model is defined as: \eqn{y_i = 1} if \eqn{z_i > 0}, and \eqn{y_i = 0} if \eqn{z_i \leq 0},
#' where \eqn{z_i = X_i \beta + \epsilon_i}, \eqn{\epsilon_i \sim N(0, \sigma^2)}.
#' A prior distribution is assumed for both \eqn{\beta} and \eqn{\sigma^2}.
#' Specifically, \eqn{\beta \sim N(\mu_{\beta}, \Sigma_{\beta})}, and \eqn{\sigma^2 \sim \text{Inv-Gamma}(\alpha_0, \beta_0)}.
#'
#' @param X A matrix of predictor variables of dimension \eqn{n \times p}, where \eqn{n} is the number of observations and \eqn{p} is the number of predictors.
#' @param y A binary vector (0/1) of response variables of length \eqn{n}, where \eqn{y_i \in \{0,1\}}.
#' @param beta A numeric vector of initial values for the regression coefficients \eqn{\beta} of length \eqn{p}.
#' @param Sigma_prior A \eqn{p \times p} matrix representing the prior covariance matrix for \eqn{\beta}.
#' @param mu_prior A numeric vector of length \eqn{p} representing the prior mean vector for \eqn{\beta}.
#' @param sigma2 A numeric value representing the initial value of \eqn{\sigma^2}.
#' @param alpha_0 A numeric value representing the shape parameter of the inverse-gamma prior for \eqn{\sigma^2}.
#' @param beta_0 A numeric value representing the scale parameter of the inverse-gamma prior for \eqn{\sigma^2}.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{beta}}{A numeric vector of updated values for the regression coefficients \eqn{\beta}.}
#'   \item{\code{sigma2}}{An updated value for the variance parameter \eqn{\sigma^2}.}
#'   \item{\code{z}}{A numeric vector of latent variables \eqn{z} sampled from the truncated normal distribution.}
#' }
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' X <- matrix(rnorm(100 * 3), 100, 3)  # 100 observations, 3 predictors
#' y <- rbinom(100, 1, 0.5)  # Binary outcome (0 or 1)
#' beta <- rep(0, 3)  # Starting beta values
#' Sigma_prior <- diag(3)  # Prior covariance matrix
#' mu_prior <- rep(0, 3)  # Prior mean for beta
#' sigma2 <- 1  # Initial value for sigma^2
#' alpha_0 <- 2  # Shape parameter for prior on sigma^2
#' beta_0 <- 2  # Scale parameter for prior on sigma^2
#'
#' result <- gibbs_binomial_probit(X, y, beta, Sigma_prior, mu_prior, sigma2, alpha_0, beta_0)
#' print(result$beta)  # Updated beta values
#' print(result$sigma2)  # Updated sigma^2 value
#'
#' @importFrom truncnorm rtruncnorm
#' @export
gibbs_binomial_probit <- function(X, y, beta, Sigma_prior,
                                  mu_prior, sigma2 = 1, alpha_0 = 2,
                                  beta_0 = 2,
                                  debug = FALSE) {

if(debug) browser()
  # Number of observations and predictors
  n <- length(y)
  p <- ncol(X)

  # Step 1: Sample latent variable z from the truncated normal distribution
  z <- rep(0, n)

  for (i in 1:n) {
    mu_i <- X[i, ] %*% beta

    if (y[i] == 1) {
      # Sample z from truncated normal distribution, truncated below at 0 for y == 1
      z[i] <- rtruncnorm(1, a = 0, mean = mu_i, sd = sqrt(sigma2))
    } else {
      # Sample z from truncated normal distribution, truncated above at 0 for y == 0
      z[i] <- rtruncnorm(1, b = 0, mean = mu_i, sd = sqrt(sigma2))
    }
  }

  # Step 2: Sample beta from the normal posterior
  Sigma_posterior <- solve(t(X) %*% X / sigma2 + solve(Sigma_prior))
  mu_posterior <- Sigma_posterior %*% (t(X) %*% z / sigma2 + solve(Sigma_prior) %*% mu_prior)

  # Draw a sample for beta from the multivariate normal distribution
  beta_new <- mvrnorm(1, mu_posterior, Sigma_posterior)

  # Step 3: Sample sigma^2 from the inverse-gamma posterior
  alpha_sigma <- alpha_0 + n / 2
  beta_sigma <- beta_0 + 0.5 * sum((z - X %*% beta_new)^2)

  # Sample sigma^2 from inverse-gamma distribution
  sigma2_new <- 1 / rgamma(1, shape = alpha_sigma, rate = beta_sigma)

  return(list(beta = beta_new, sigma2 = sigma2_new, z = z))
}
