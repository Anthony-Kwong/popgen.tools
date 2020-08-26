#' pseudo_hap function
#' 
#' Simulates pseudo haplotypes for a genome matrix. The matrix is broken into pairs of 
#' rows. For each row pair, we randomly sample an element for each column to make
#' new pseudo haplotypes. Function outputs a matrix of pseudo haplotypes with half
#' the number of rows as the original genome matrix.
#'
#' @param G : A binary genome matrix: The number of rows must be even. 
#' @param seed: A numeric random seed.
#'
#' @return : A smaller genome matrix consisting of the pseudo-haplotypes of G. Output has
#' half the number of rows as G. 
#' @export
#'
#' @examples
#' 
#' @importFrom tester is_not_even
#' @importFrom rlist list.rbind
pseudo_hap <- function (G, seed = NA){

  nsam = nrow(G)
  if(tester::is_not_even(nsam)){
    stop("Number of rows in G must be even. They are odd right now.")
  }
  if(is.na(seed)){
    seed = sample.int(.Machine$integer.max, size = 1)
  }
  set.seed(seed)
  p_seeds = sample.int(.Machine$integer.max, size = nsam/2)
  
  snp = ncol(G)
  nsam = nrow(G)
  #find index of the first row for each haplotype pair
  hap_indices = seq(1, nsam, by = 2)
  pseu_hap = list()
  length(pseu_hap) = nsam/2
  
  for(k in 1:length(pseu_hap)){
    set.seed(p_seeds[k])
    choice = sample(c(0,1) , replace=TRUE, size= snp)
    new_hap = rep(NA, snp)
    for(i in 1:snp){
      if(choice[i] == 1){
        new_hap[i] = G[ hap_indices[k], i]
      } else {
        new_hap[i] = G[hap_indices[k], i]
      }
    }
    pseu_hap[[k]] = new_hap
  }
  #bind all the rows together into a matrix
  final = rlist::list.rbind(pseu_hap)
  
  return(final)
}
