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

#' Parse a complex formula with optional '|' symbol
#'
#' This function parses a formula that may contain a '|' symbol, indicating two
#' separate formulas. It handles cases where there is no '|' symbol, in which case
#' it returns the original formula. If more than one '|' symbol is found, the
#' function throws an error.
#'
#' @param formula A formula object or a string that can be converted to a formula.
#' It can optionally contain a '|' symbol separating two formulas.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{left_formula}{The formula to the left of the '|' or the original formula if no '|' is present.}
#'   \item{right_formula}{The formula to the right of the '|' or NULL if no '|' is present.}
#' }
#'
#' @examples
#' # Example 1: Formula with one '|'
#' f.X <- "y ~ x + b | x ~ g + v"
#' parsed_formulas <- parse_complex_formula(f.X)
#' print(parsed_formulas$left_formula)  # Should print: y ~ x + b
#' print(parsed_formulas$right_formula)  # Should print: x ~ g + v
#'
#' # Example 2: Formula without '|'
#' f.X_no_pipe <- "y ~ x + b"
#' parsed_formulas_no_pipe <- parse_complex_formula(f.X_no_pipe)
#' print(parsed_formulas_no_pipe$left_formula)  # Should print: y ~ x + b
#' print(parsed_formulas_no_pipe$right_formula)  # Should print: NULL
#'
#' # Example 3: Formula with more than one '|'
#' f.X_multiple_pipes <- "y ~ x + b | x ~ g + v | z ~ w"
#' # This will throw an error
#' # parsed_formulas_multiple_pipes <- parse_complex_formula(f.X_multiple_pipes)
#'
#' @export
parse_complex_formula <- function(formula) {
  # Check if the input is a character string
  if (is.character(formula)) {
    # Try to validate and convert the string to a formula
    tryCatch({
      formula <- as.formula(formula)
    }, error = function(e) {
      stop("Error: Invalid formula string. Could not parse.")
    })
  }
  # Deparse the formula object to a string for processing
  formula_str <- paste(deparse(formula, width.cutoff = 500), collapse = "")

  # Use gregexpr to find the positions of '|' symbols
  pipe_positions <- gregexpr("\\|", formula_str)[[1]]

  # Count the occurrences of '|' (ignore if the result is -1)
  pipe_count <- sum(pipe_positions != -1)

  # If more than one '|' is found, throw an error
  if (pipe_count > 1) {
    stop("Error: The formula contains more than one '|' symbol.")
  }

  # Check if the formula contains exactly one '|'
  if (pipe_count == 1) {
    # Split the formula into two parts using the '|' as the separator
    formula_parts <- strsplit(formula_str, "\\|")[[1]]

    # Parse both parts as complete formulas
    left_formula <- as.formula(formula_parts[1])  # Left part of the '|'
    right_formula <- as.formula(formula_parts[2]) # Right part of the '|'

    # Return both formulas as a list
    return(list(left_formula = left_formula, right_formula = right_formula))
  } else {
    # If there's no '|', return the original formula in a list with NULL for right part
    return(list(left_formula = formula, right_formula = NULL))
  }
}



