#' add_missingness function
#' 
#' Randomly adds missingness to a genome matrix. Each element in the matrix becomes NA with
#' a fixed probability (missing_rate). The trials are independent for each matrix element. 
#' 
#'
#' @param G : A binary, numeric genome matrix. 
#' @param missing_rate : Probability for any given element of G to become NA.
#' @param seed : Optional. A random seed for reproducibility. 
#'
#' @return A genome matrix with NA elements.
#' @export
#'
#' @examples   G = matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
#' add_missingness(G)
#' @importFrom purrr rbernoulli is_empty


add_missingness <- function(G, missing_rate, seed = NA){
  
  #check input arguments
  if(missing_rate > 1 || missing_rate < 0){
    stop("missing rate must be a numeric between 0 and 1")
  }
  if(length(missing_rate) != 1){
    stop("missing rate must be a single number, not a vector.")
  }
  
  if(is.na(seed)){
    seed = sample.int(.Machine$integer.max, size = 1) 
  }
  
  #computing basic variables
  snp = ncol(G)
  nsam = nrow(G)
  
  #setting random seeds for reproducibility
  set.seed(seed)
  NA_seeds = sample.int(.Machine$integer.max, size = nsam) 
  
  for(i in 1:nsam){
    set.seed(NA_seeds[i])
    index = purrr::rbernoulli(n=snp, p = missing_rate)
    G[i,which(index == TRUE)] = NA
  }
  
  return(G)
}
