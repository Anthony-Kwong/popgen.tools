#' sample_vec function
#' 
#' Randomly samples an element from a vector. If the vector contains only one element, that
#' element is returned. 
#'
#' @param x: A vector.
#'
#' @return An random element selected from the input vector.
#' @export
#'
#' @examples x <- c(1,3,4,5)
#' sample_vec(x)
sample_vec <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}