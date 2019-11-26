#' norm_vec function
#' 
#' Normalizes a numeric vector by dividing it by the sum of all its elements.
#'
#' @param vec: a numeric vector
#'
#' @return a normalised numeric vector
#' @export
#'
#' @examples x<-rep(1:10)
#' norm_vec(x)
norm_vec<-function(vec){
  return(vec/sum(vec))
}