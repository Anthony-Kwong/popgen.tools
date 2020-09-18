#' zero_impute function
#' 
#' Changes all NAs in a matrix into 0. A conservative approach as it assumes that 
#' there is no derived allele present where there is missing information.
#'
#' @param G: A NumericMatrix with NAs. 
#'
#' @return A new matrix where all the NA elements of G have been replaced by 0.
#' @export
#'
#' @examples G <- rbind(c(1, 0, 1, NA, 0), c(1, NA, NA, 0, 1), c(0, 0, 1, 0, NA))
#'  zero_impute(G)
zero_impute <- function (G){
  G[is.na(G)] = 0
  return(G)
}