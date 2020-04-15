#' is_genome_matrix function
#' 
#' Checks if an R object is a valid genome matrix. A genome matrix is a binary 
#' numeric matrix consisting of only 0's and 1's. Each column must have at least
#' one entry which is 1. 
#'
#' @param x: an R object 
#'
#' @return a logical indicating whether the input is a valid genome matrix
#' @export
#'
#' @examples seq <-matrix(sample(0:1, size =16 , replace = TRUE), nc = 4)
#' is_genome_matrix(seq)
#' @importFrom tester is_numeric_matrix
is_genome_matrix <- function (x){
  #allow null windows to pass with warning
  if(any(is.na(x))){
    warning("NAs found in block. Probably null window. Summary stats will give NAs for this block")
    return(T)
  }
  
  
  if(tester::is_numeric_matrix(x) == F){
    return (F)
  } else if (all (x %in% c(0, 1) ) ==F){
    return (F)
  } else if (any ( apply(x, 2, function(vec) { all(vec == vec[1]) }) ) ){
    #this condition checks if any column has all the same values
    return (F)
  }  else {
    return (T)
  }
  
}
