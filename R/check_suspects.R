#' Prepares data
#'
#' Function to prepare data for an object to be taken in to Gibbs sampler
#'
#' @param dir A directory containing results
#' @return A vector of file names that are suspect
#'
#' @export
check_suspects <- function(dir){

  # find the files
  flist <- list.files(path = dir,
                      pattern = "rds")
  suspect <- NULL
  for(fh in flist){
    .a <- readRDS(paste0(dir,fh))
    w <- .a$res$weights
    # check
    if(max(colSums(w))>1.1) {
      suspect <- rbind(suspect,c(fh,
                                 round(min(colSums(w)),4),
                                 round(max(colSums(w)),4)))
      } else {
        if(max(w) < 0.0001) suspect <- rbind(suspect,c(fh,round(min(w),4),round(max(w),4)))
        if(max(w) > 10) suspect <- rbind(suspect,c(fh,round(min(w),4),round(max(w),4)))
    }
  }

  return(suspect)
}
