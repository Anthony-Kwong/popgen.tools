#' fixed_clus_hap function
#' 
#' Takes a genome matrix and clusters the rows into n_clus clusters. If the 
#' number of unique rows is higher than n_clus, then the number of unique rows
#' is used instead. 
#'
#' @param G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual. 
#' @param n_clus: Number of clusters to cluster the rows of G. 
#'
#' @return A Numeric vector indicating the clusters to which each row belongs
#' @export
#'
#' @examples G = matrix(sample(0:1, size = 40 , replace = TRUE), nc = 8)
#' fixed_clus_hap(G, 3)
#' @importFrom tester is_integer

fixed_clus_hap <- function (G, n_clus){
  #check inputs
  
  if(tester::is_integer(n_clus)==F){
    stop("max_clus must be an integer.")
  }
  if(n_clus > nrow(G)){
    stop("max_clus must not exceed the number of rows in G")
  }
  
  #find the number of unique rows in genome matrix. 
  #max_clus must not exceed the number of unique rows as this would break clustering
  unique_dp = sum(!duplicated(G))
  
  if(unique_dp < n_clus){
    n_clus = unique_dp
  }
  
  #cluster rows
  G_clus = kmeans(G, centers = n_clus)
  clus_vec = G_clus$cluster
  
  return(clus_vec)
}