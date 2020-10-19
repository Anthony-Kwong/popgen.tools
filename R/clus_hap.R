#' clus_hap function
#' 
#' Takes a genome matrix and clusters the rows. The best number of clusters is chosen
#' by having the highest silhouette value. Returns a clustering vector. The i^th element
#' of the cluster vector is the cluster to which row i has been assigned. 
#'
#'
#' @param G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.
#' @param max_clus : Integer for the maximum number of clusters to consider.
#'
#' @return A Numeric vector indicating the clusters to which each row belongs
#' @examples G = matrix(sample(0:1, size =25 , replace = TRUE), nc = 5)
#' clus_hap(G,3)
#' @export
#' @importFrom tester is_integer
#' @importFrom factoextra fviz_nbclust
clus_hap <- function(G, max_clus){
  #check inputs
  if(is_genome_matrix(G)==F){
    stop("invalid genome matrix.")
  }
  if(tester::is_integer(max_clus)==F){
    stop("max_clus must be an integer.")
  }
  if(max_clus > nrow(G)){
    stop("max_clus must not exceed the number of rows in G")
  }
  
  #tune the number of clusters
  tune_clus = factoextra::fviz_nbclust(G, kmeans, method = "silhouette",
                                       k.max = max_clus)
  sil_data = tune_clus$data
  #take the number of clusters with the highest silhouette value
  best_clus = as.numeric(sil_data$clusters[which.max(sil_data$y)])
  
  #cluster rows
  G_clus = kmeans(G, centers = best_clus)
  clus_vec = G_clus$cluster
  
  return(clus_vec)
}