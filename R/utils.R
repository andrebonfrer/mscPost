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

#' Parse a complex formula with optional '|' symbols
#'
#' This function parses a formula that may contain one or more '|' symbols, indicating
#' multiple formulas. It returns a list of formula objects split by the '|' symbol.
#' If no '|' symbol is found, it returns the original formula as a single list element.
#'
#' @param formula A formula object or a string that can be converted to a formula.
#' It can optionally contain one or more '|' symbols separating formulas.
#'
#' @return A list of formula objects.
#'
#' @examples
#' # Example 1: Formula with multiple '|'
#' f.X <- "y ~ x + b | x ~ g + v | z ~ w"
#' parsed_formulas <- parse_complex_formula(f.X)
#' print(parsed_formulas)  # Should print the list of formulas: y ~ x + b, x ~ g + v, z ~ w
#'
#' # Example 2: Formula without '|'
#' f.X_no_pipe <- "y ~ x + b"
#' parsed_formulas_no_pipe <- parse_complex_formula(f.X_no_pipe)
#' print(parsed_formulas_no_pipe)  # Should print: y ~ x + b
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

  # Split the formula string by the '|' symbol
  formula_parts <- strsplit(formula_str, "\\|")[[1]]

  # Trim whitespace from each part and convert each to a formula
  formula_list <- lapply(formula_parts, function(part) {
    trimmed_part <- trimws(part)
    as.formula(trimmed_part)
  })

  # Return the list of formulas
  return(formula_list)
}



