#' no_scientific function
#'
#'Changes a number away from scientific notation. Discoal doesn't like scientific notation. 
#'
#' @param x some numeric scalar
#'
#' @return the number x as a character vector without scientific notation
#' @export
#'
#' @examples
#' change 1e6 to 1000000
#' num<-no_scientific(1e6)
#'
no_scientific<-function(x){
  num=format(x,scientific = F, trim = T)
  return (num)
}
