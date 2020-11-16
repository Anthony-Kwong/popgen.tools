#' clus_hap function
#' 
#' Takes a genome matrix and clusters the rows. The best number of clusters (>=3)
#' is chosen by having the highest silhouette value. Returns a clustering vector. The i^th 
#' element of the cluster vector is the cluster to which row i has been assigned. 
#'
#'
#' @param G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.
#' @param max_clus : Integer for the maximum number of clusters to consider. If the number of 
#' unique rows is smaller than max_clus, the number of unique rows will be used instead. Must 
#' be at least 3.
#'
#' @return A Numeric vector indicating the clusters to which each row belongs
#' @examples G = matrix(sample(0:1, size =25 , replace = TRUE), nc = 5)
#' clus_hap(G,3)
#' @export
#' @importFrom tester is_integer
#' @importFrom factoextra fviz_nbclust
clus_hap <- function(G, max_clus){
  #check inputs
  # if(is_genome_matrix(G)==F){
  #   stop("invalid genome matrix.")
  # }
  
  if(tester::is_integer(max_clus)==F){
    stop("max_clus must be an integer.")
  }
  if(max_clus > nrow(G)){
    stop("max_clus must not exceed the number of rows in G")
  }
  if(max_clus < 3){
    stop("max_clus must be at least 3 to enable computation of h123.")
  }
  
  #find the number of unique rows in genome matrix. 
  #max_clus must not exceed the number of unique rows as this would break clustering
  unique_dp = sum(!duplicated(G))
  
  if(unique_dp < max_clus){
    max_clus = unique_dp
  }
  
  #tune the number of clusters
  tune_clus = factoextra::fviz_nbclust(G, kmeans, method = "silhouette",
                                       k.max = max_clus)
  sil_data = tune_clus$data
  #take the number of clusters with the highest silhouette value
  best_clus = as.numeric(sil_data$clusters[which.max(sil_data$y)])
  
  #ensure we have at least 3 clusters so we don't get NAs for h123
  if(best_clus < 3){
    best_clus = 3
  }
  
  #cluster rows
  G_clus = kmeans(G, centers = best_clus)
  clus_vec = G_clus$cluster
  
  return(clus_vec)
}