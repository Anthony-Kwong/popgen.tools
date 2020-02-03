#' norm_vec function
#' 
#' Normalizes a numeric vector by taking away the mean and dividing by the standard deviation, for each element. 
#'
#' @param vec: a numeric vector
#' @return a normalised numeric vector
#' @export
#' @examples x<-rep(1:10)
#' norm_vec(x)
norm_vec<-function(vec){
  e<-mean(vec)
  std<-var(vec)^0.5
  return( (vec-e)/std )
}