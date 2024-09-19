#' Prepares data
#'
#' Function to prepare data for an object to be taken in to Gibbs sampler
#'
#' @param dir A directory containing results
#' @return A vector of file names that are suspect
#'
#' @export
check_suspects <- function(dir){

  # point to the right directory
  setwd(dir)

  flist <- list.files(path = dir,
                      pattern = "rds")
  suspect <- NULL
  for(fh in flist){
    .a <- readRDS(paste0(dir,fh))
    w <- .a$res$weights
    # check
    if(max(w) < 0.01) suspect <- rbind(suspect,c(fh,min(.a),max(.a)))
    if(max(w) > 10) suspect <- rbind(suspect,c(fh,min(w),max(w)))
  }

  return(suspect)
}
