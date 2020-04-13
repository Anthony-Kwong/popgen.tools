#' downsample_mat function
#' 
#' Downsample a NumericMatrix by randomly sampling the columns without replacement. 
#'
#' @param G : A NumericMatrix designating a binary genome matrix consisting of 1's and 0's. 
#' @param p : Proportion of columns to downsample. 
#' @param seed: Random seed for sampling columns. (optional)
#'
#' @return A smaller NumericMatrix with randomly selected columns from G.
#' @export
#'
#' @examples seq <-matrix(sample(0:1, size =40 , replace = TRUE), nc = 10)
#' downsample_mat(seq, 0.2)
#' @importFrom tester has_NA
downsample_mat = function (G , p, seed = NA){
  
  if(tester::has_NA(G)){
    warning("Input matrix is NA")
    return(NA)
  }
  
  if(is.matrix(G)!=T){
    stop("Parameter G must be a matrix.")
  }
  
  if(p>1 || p<=0){
    stop("Parameter p must be between 0 and 1. It is a proportion.")
  }
  
  cols = ncol(G)
  n = round( ncol(G)*p )  
  
  if(is.na(seed)){
    seed = sample.int(.Machine$integer.max,1)
  }
  
  set.seed(seed)
  sam_col = sample(1:cols, n)
  indices = sort(sam_col)
  
  G2 = G[,indices]
  
  return(G2)
}

