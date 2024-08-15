#' Winsorizing function
#'
#' This function runs a Gibbs sampler for a given set of parameters.
#'
#' @param x A variable name.
#' @param probs A numeric list of lower and upper tails
#' @export

# utilities
winsorize <- function(x, probs = c(0.025, 0.975)) {
  # Ensure the probs are in the correct range
  if (any(probs < 0 | probs > 1)) {
    stop("Probabilities must be between 0 and 1")
  }

  # Compute the quantiles
  q <- quantile(x, probs = probs, na.rm = TRUE)

  # Winsorize the values
  x_winsorized <- pmax(pmin(x, q[2]), q[1])

  return(x_winsorized)
}

#' Efficient draw from multivariate normal
#'
#' Useful for very large draws of a multivariate vector
#' uses a Cholesky decomposition
#'
#' @param n Number of draws.
#' @param mu A vector of means.
#' @param Sigma A covariance matrix
#' @importFrom Matrix Matrix chol
#' @export
rMVNormCovariance <- function(n, mu, Sigma){
  p <- length(mu)
  Z <- Matrix::Matrix(rnorm(p*n), nrow = p, ncol = n)
  L <- Matrix::t(Matrix::chol(Sigma)) # By default R's chol function returns upper cholesky factor
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
  return(X)
}

