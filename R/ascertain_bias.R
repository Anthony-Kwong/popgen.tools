#' ascertain_bias function
#' 
#' Function to simulate ascertainment bias in ancient DNA. Takes in a genome matrix G
#' where some rows represent sequences from an outgroup (usually 2). The function 
#' identifies all heterozygous sites in the outgroup, then subsets the G at these 
#' columns. The outgroup is excluded in the output matrix. 
#'
#' @param G : A binary genome matrix. 
#' @param index : A vector of indices indicating the rows of G which represent the outgroup.
#'
#' @return A genome matrix consisting of the columns of G where the outgroup was hetereozygous.
#' The rows belonging to the outgroup are also excluded. 
#' @export
#'
#' @examples G <-matrix(sample(0:1, size = 40, replace = TRUE), nc = 8)
#' ascertain_bias(G, c(7,8))
ascertain_bias <-function(G, index){
  outgroup_G = G[index,]
  het_sites = het_finder(outgroup_G)
  final_G = matrix_subset(G,het_sites)[-index,]
  return(final_G)
}