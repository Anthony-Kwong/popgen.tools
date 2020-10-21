#' majority_flip function
#' 
#' Takes a binary NumericMatrix of 0's and 1's. For each column, one entry from the
#' minority group is randomly selected to be switched to the majority.
#'
#' @param G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual. 
#' @param seed: A random seed.
#' @return A NumericMatrix where one element has been randomly switched in each column.
#' @export
#'
#' @examples G = matrix(sample(0:1, size =20 , replace = TRUE), nc = 5)
#' majority_flip(G)
majority_flip <- function (G, seed){
  #check inputs
  
  
  snp = ncol(G)
  #generate random seeds
  seed_vec = sample.int(.Machine$integer.max, size = snp)
  
  #loop over each column
  for(c in 1:snp){
    sel_col = G[,c]
    num_0 = count(sel_col, 0)
    num_1 = count(sel_col, 1)
    
    if(num_0 == num_1){
      next
    } else if (num_0 > num_1){
      pos_flip = which(sel_col == 1)
      set.seed(seed_vec[c])
      flip = sample_vec(pos_flip)
      #print(flip)
      G[flip,c] = 0
    } else {
      pos_flip = which(sel_col == 0)
      set.seed(seed_vec[c])
      flip = sample_vec(pos_flip)
      #print(flip)
      G[flip,c] = 1
    }
  }
  return(G)
}