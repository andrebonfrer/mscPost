#' Prepares data
#'
#' Function to prepare data for an object to be taken in to Gibbs sampler
#'
#' @param dta A data object containing outcomes and covariates
#' @param res An object for modified augmented synthetic control results
#' produced by augsynth multisynth
#' @param f.X A formula for X - must include "tvg.dummy" for intervention
#' @param f.Z A formula for moderators for intervention HTE
#' @param flags A set of needed parameters
#' @param verbose Whether to print a lot of information as it's running
#' @importFrom stats as.formula binomial dnorm glm model.matrix pnorm predict quantile rgamma rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Matrix Diagonal Matrix
#' @importFrom data.table := is.data.table uniqueN as.data.table
#'
# makes X, y, Z, and W as a list
#' @export
prepare_data <- function(dta = NULL, res = NULL,
                         f.X = as.formula(),
                         f.Z = as.formula(),
                         flags = NULL,
                         verbose = TRUE){

  utils::globalVariables(c("tvg.dummy", "id", "V1", "GR", "SCweight", "..xvars", "dv", "wID"))

  # check for dataset and weights
  if(is.null(dta) || !is.data.table(dta)) stop("Provide a data table!")
  if(is.null(res)) stop("Provide a results matrix!")

  # Define spinner symbols
  spinner <- c("|", "/", "-", "\\")

  f.X <- as.formula(f.X)
  f.Z <- as.formula(f.Z)

  # parse the y~X formula
  dvname <- all.vars(f.X)[1]

  # selection equation code below - later make as (part of) formula
  # first stage object for selection equation
  if(verbose) cat("Running first stage (selection) equation. Be patient.\n")
  f1 <- paste0("tvg.dummy ~ factor(wID) + ",
               paste0(flags$columns_instruments, collapse = "+"),
               "+",
               paste0(flags$cov.var.first.stage, collapse="+"),
               "+ email_count_week + I(email_count_week^2) + count.push"
  )

  # Fit the selection equation (first stage)
  first_stage <- glm(f1 ,
                     family = binomial(link = "probit"),
                     data = dta)

  # add the generalized residual
  PR <- predict(first_stage, type = "response")
  pdf_values <- dnorm(PR)
  cdf_values <- pnorm(PR)
  # Compute the Inverse Mills Ratio
  imr <- pdf_values / (1 - cdf_values)

  # Add Generalized Residual to data set
  dta$GR <- dta$tvg.dummy*imr + (1-dta$tvg.dummy)*imr

  # finished with selection equation set up

  # weights matrix from data set
  weightmatrix <- res$weights

  # Find all N customers
  idlist <- unique(dta$id)

  # this part figures out how many treated there are
  # just those treated (goal setters in obs period)
  tlist <- dta[,max(tvg.dummy),by=id][V1==1, id]

  J <- length(idlist)   # Number of individuals including never treated
  J0 <- length(tlist)
  if(verbose) cat("A total of", J, "customers.\n")
  if(verbose) cat("There are:", J0, "treated individuals.\n")

  # Define parameters
  T <- uniqueN(dta$wID)
  N <- (1+J-J0)*T  # Number of observations per outcome

  # Create predictor matrices X
  if(verbose) cat("Creating y, X, and w matrix. (This takes a while.)\n")
  y_list <- X_df_list <- X_list <- weights_list <- list()
  for(j0 in 1:J0) {
    cat("\r", spinner[(j0 %% length(spinner)) + 1], sep = "")
    wj <- weightmatrix[,j0]
    # find the treated individual in the donor pool of all N
    # needed? I think we keep it at zero?
    j.idx <- match(tlist[j0], idlist)
    wj[j.idx] <- 1

    # those that are in donor pool for this person, including this person
    .idlist <- c(tlist[j0], setdiff(idlist,tlist))

    # cv for X (set up outside loop)
    cXv <- model.matrix(f.X,data=dta)
    # set up panel for one person versus donor, aliging weights
    .a <- as.data.table(cbind(cXv,dta[,list(id,dv = get(dvname),GR)]))
    .a[,SCweight := wj[match(id,idlist)]]

    # check
    .sum <- .a[SCweight>0, sum(tvg.dummy)]
    if(.sum==0) stop(paste0("Something is wrong, sum of tvg.dummy is zero for index: ", j0))

    # get the X covariates for each individual
    xvars <- setdiff(colnames(.a), c("id", "SCweight","dv"))
    .aa <- .a[id %in% .idlist,..xvars]

    X_list[[j0]] <- as.matrix(.aa)
    X_df_list[[j0]] <- cbind(1, .aa[,list(tvg.dummy,GR)])

    y_list[[j0]] <- .a[id %in% .idlist,dv]

    # Generate weights (for simplicity, we use random positive values)
    weights_list[[j0]] <- .a[id %in% .idlist,SCweight]

  }

  # figure out K from here
  K <- ncol(X_list[[1]])

  # Combine X_list, y_list, and weights_list into arrays
  X_array <- array(NA, dim = c(J0, N, K))
  y_array <- array(NA, dim = c(J0, N))
  y_long <- weights_array <- rep(0, length=J0*N)

  for (j0 in 1:J0) {
    X_array[j0, , ] <- X_list[[j0]]
    y_array[j0, ] <- y_list[[j0]]
    weights_array[(j0-1)*N+1:N] <- weights_list[[j0]]
    y_long[(j0-1)*N+1:N] <- y_list[[j0]]
  }

  Xb <- Matrix::bdiag(X_list)

  W <- Diagonal(x=weights_array)

  Y <- Matrix(data=y_long, byrow=TRUE)

  Z.dm <- model.matrix(f.Z, data = dta[wID==T & id %in% tlist,])

  .cn <- colnames(.aa)

  cat("Finished with creation of data object.\n")

  return(list(Y_block = Y,
              X_block = Xb,
              W = W,
              Z = Z.dm,
              cov = list(Xcols = .cn,
                         Zcols = colnames(Z.dm),
                         intX = which(.cn=="tvg.dummy"))
  ))
}
