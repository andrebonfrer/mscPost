# specify any global variables here
utils::globalVariables(c("tvg.dummy", "id", "V1", "GR", "SCweight", "..xvars", "dv", "wID"))

#' Prepares data
#'
#' Function to prepare data for an object to be taken in to Gibbs sampler
#'
#' @param dta A data object containing outcomes and covariates
#' @param res An object for modified augmented synthetic control results
#' produced by augsynth multisynth
#' @param f.X A formula for X - must include "tvg.dummy" for intervention
#'
#' @param f.Z A formula for moderators for intervention HTE. For example,
#' we have a formula f.Z such as:
#' tvg.dummy ~ 1 + Z1 + Z2 + Z3 | Z1 ~ Q1 + Q2 | Z2 ~ Q3 + Q4
#' and this will have instruments Q1 and Q2 for Z1, and Q3 and Q4 for Z2.
#'
#' @param flags A set of needed parameters
#' @param verbose Whether to print a lot of information as it's running
#' @importFrom stats as.formula binomial dnorm glm model.matrix pnorm predict quantile rgamma rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Matrix Diagonal Matrix
#' @importFrom data.table := is.data.table uniqueN as.data.table
#' @export
prepare_data <- function(dta = NULL, res = NULL,
                         f.X = as.formula(),
                         f.Z = as.formula(),
                         flags = NULL,
                         verbose = TRUE){


  # check for dataset and weights
  if(is.null(dta) || !is.data.table(dta)) stop("Provide a data table!")
  if(is.null(res)) stop("Provide a results matrix!")

  # Define spinner symbols
  spinner <- c("|", "/", "-", "\\")

  # first stage parts
  f.X.parsed <- parse_complex_formula(f.X)

  # second stage equation
  f.Z.parsed <- parse_complex_formula(f.Z)

  # parse the y~X formula (LHS of the first formula object for f.X)
  dvname <- all.vars(f.X.parsed[[1]])[1]

  # selection equation code below - later make as (part of) formula
  # first stage object for selection equation
  if(!is.null(f.X.parsed[[2]])) {
    if(verbose) cat("Running first stage (selection) equation. Be patient.\n")

    # Fit the selection equation (first stage)
    first_stage <- glm(f.X.parsed[[2]],
                     family = binomial(link = "probit"),
                     data = dta)

    # add the generalized residual
    PR <- predict(first_stage, type = "response")

    # Calculate and add Generalized Residual (GR) to data set
    dta$GR <- dta$tvg.dummy*dnorm(PR)/pnorm(PR) -
                (1-dta$tvg.dummy)*dnorm(-PR)/pnorm(-PR)

    # save the model matrix here
    GRX <- model.matrix(f.X.parsed[[2]], data = dta)
    GRY <- dta$tvg.dummy
  }  # finished with selection equation set up (if included)

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
  X_idlist <- NULL
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
    cXv <- model.matrix(f.X.parsed[[1]], data=dta)
    # set up panel for one person versus donor, aligning weights
    if(!is.null(f.X.parsed[[2]])) {
      .a <- as.data.table(cbind(cXv,dta[,list(id,dv = get(dvname),GR)]))
    } else {
      .a <- as.data.table(cbind(cXv,dta[,list(id,dv = get(dvname))]))
    }
    .a[,SCweight := wj[match(id,idlist)]]

    # check
    .sum <- .a[SCweight>0, sum(tvg.dummy)]
    if(.sum==0) stop(paste0("Something is wrong, sum of tvg.dummy is zero for index: ", j0))

    # get the X covariates for each individual
    xvars <- setdiff(colnames(.a), c("id", "SCweight","dv"))
    .aa <- .a[id %in% .idlist,..xvars]

    X_list[[j0]] <- as.matrix(.aa)
    if(!is.null(f.X.parsed[[2]])) {
      X_df_list[[j0]] <- cbind(1, .aa[,list(tvg.dummy,GR)])
    } else {
      X_df_list[[j0]] <- cbind(1, .aa[,list(tvg.dummy)])
    }

    # store wID and id for matching this within sampler
    X_idlist <- rbind(X_idlist, dta[id %in% .idlist,list(id,wID)])

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

  .Zdt <- dta[wID==T & id %in% tlist,]
  Z.dm <- model.matrix(f.Z.parsed[[1]], data = .Zdt)

  # check for instruments - return NULL if it does not exist
  if(length(f.Z.parsed)>1) {
    Z.im <- list()
    dv <- NULL
      for(l in 2:length(f.Z.parsed)) {
        dv <- c(dv, all.vars(f.Z.parsed[[l]])[1])
        Z.im[[l-1]] <- model.matrix(f.Z.parsed[[l]], data = .Zdt)
      }
    Z.instruments <- list(dv = dv,
                          Z.im = Z.im)
    } else Z.instruments <- NULL

  .cn <- colnames(.aa)

  cat("Finished with creation of data object.\n")

  return(list(Y_block = Y,
              X_block = Xb,
              W = W,
              Z_block = Z.dm,
              Z.instruments = Z.instruments,
              GRX = GRX,
              GRY = GRY,
              X_idlist = X_idlist,
              dtaidx = dta[,list(id,wID)],
              cov = list(Xcols = .cn,
                         Zcols = colnames(Z.dm),
                         intX = which(.cn == all.vars(f.Z.parsed[[1]])[1]),
                         J0 = J0, J = J, T = T)
  ))
}
