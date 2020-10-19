#' clus_hstats function
#' 
#' Takes in a cluster vector of a genome matrix and computes the haplotype statistics.
#'
#' @param clus_vec: A numeric vector indicating to which cluster, each row in the original 
#' genome matrix was assigned. The clusters are labelled from 1,2,3.... up to the some number 
#' of clusters.
#'
#' @return A numeric vector of haplotype statistics (h1,h2,h12,h123).
#' @export
#'
#' @examples G = matrix(sample(0:1, size =25 , replace = TRUE), nc = 5)
#' x = clus_hap(G,3)
#' clus_hstats(x)
clus_hstats <- function (clus_vec){
  num_hap = length(unique(clus_vec))
  #compute proportions of each cluster
  hap_counts = sapply( seq(1:num_hap) , function(d){count(clus_vec,d)} )
  hap_freq = hap_counts/length(clus_vec)
  
  h1 = sum(hap_freq*hap_freq)
  h2 = h1 - hap_freq[1]^2
  h12 = h1 + 2*hap_freq[1]*hap_freq[2]
  #if we don't have at least 3 clusters. This value becomes NA.
  h123 = h12 + 2*hap_freq[1]*hap_freq[3] + 2*hap_freq[2]*hap_freq[3] 
  
  hstat = c(h1,h2,h12,h123)
  return(hstat)
}