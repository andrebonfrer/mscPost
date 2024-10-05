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


#' Replace k-th column in each block of a block diagonal matrix with CR
#'
#' This function replaces the k-th column in each block of a block diagonal matrix `M`
#' with corresponding sub-vectors from a column vector `CR`. The matrix `M` is assumed
#' to be composed of `N` blocks, each of size `T x K`. The vector `CR` is of length
#' `N * T`, and is divided into `N` sub-vectors of length `T`, each of which replaces
#' the k-th column of one block in `M`.
#'
#' @param M A block diagonal matrix with `N` blocks, each of size `T x K`. The matrix is expected
#' to be in dense format.
#' @param GRX A numeric column vector of length `N * T` to replace the k-th column in each block of the matrix.
#' @param T An integer specifying the number of rows in each block of the matrix.
#' @param K An integer specifying the number of columns in each block of the matrix.
#' @param N An integer specifying the number of blocks in the matrix.
#' @param k An integer specifying the index of the column (in each block) to be replaced with the corresponding sub-vector from `GRX`.
#' This value must be in the range \eqn{1 \leq k \leq K}.
#'
#' @return A matrix `M` with the k-th column of each block replaced by the respective sub-vectors from `GRX`.
#'
#' @examples
#' # Example of using the replace_kth_column_block_diagonal function
#' T <- 5  # Number of rows in each block
#' K <- 3  # Number of columns in each block
#' N <- 4  # Number of blocks
#'
#' # Create a block diagonal matrix M with N blocks, each of size T x K
#' M <- Matrix::bdiag(replicate(N, matrix(rnorm(T * K), T, K), simplify = FALSE))
#' M <- as.matrix(M)  # Convert the block diagonal matrix to a dense format
#'
#' # Column vector CR to replace the k-th column in each block (length N * T)
#' CR <- rnorm(N * T)
#'
#' # Replace the 2nd column (k = 2) of each block with the corresponding sub-vectors from GRX
#' k <- 2
#' M_updated <- replace_kth_column_block_diagonal(M, CR, T, K, N, k)
#'
#' print(M_updated)  # Print the updated matrix
#'
#' @export
replace_kth_column_block_diagonal_fast <- function(M, CR, T, K, N, k) {
  # Check that CR has the correct length
  if (length(CR) != N * T) {
    stop("CR must be a column vector of length N * T.")
  }

  if (k < 1 || k > K) {
    stop("k must be a valid column index between 1 and K.")
  }

  # Find the indices of the k-th column in each block
  for (i in 1:N) {
    row_start <- (i - 1) * T + 1
    row_end <- i * T
    col_index <- (i - 1) * K + k  # The column index in the large matrix M

    # Replace the values in the k-th column of the i-th block with the sub-vector from GRX
    M[row_start:row_end, col_index] <- CR[row_start:row_end]
  }

  return(M)
}

