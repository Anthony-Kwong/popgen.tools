#' age_DNA function
#' 
#' A function to simulating DNA aging on discoal simulations. To simulate missingness, for each row 
#' a percentage of elements are randomly sampled and become missing. To simulate deamination
#' a percentage of columns are randomly sampled according to trans_rate, to become transition sites. 
#' For each transition site we flip a fair coin. If heads we change 0's to 1's, with probability as 
#' specified by the dmg_rate. If tails, we change 1's to 0's, with probability as 
#' specified by the dmg_rate.
#' 
#' This method is adapted from doi: 10.1534/genetics.112.139949
#'
#' @param G : A binary, numeric genome matrix from discoal. 
#' @param missing_rate : Probability of elements that are randomly sampled to become NA for each row in G.
#' @param trans_prop: Proportion of columns (i.e. sites) that are chosen to be transition sites.
#' Default is 0.776.
#' @param dmg_rate: Probability of a element in a transition column changing from 0 to 1, or 1 to 0. 
#' Default is 0.05. 
#' @param seed: A numeric random seed. Optional.
#' @return
#' @export
#'
#' @examples
#' 
#' @importFrom purrr rbernoulli is_empty
age_DNA <- function(G, missing_rate, trans_rate = 0.776, dmg_rate = 0.05, seed = NA){

  #check inputs
  if(is_genome_matrix(G) == F){
    stop("G is not a valid genome matrix.")
  }
  if(missing_rate > 1 || missing_rate < 0){
    stop("missing rate must be a numeric between 0 and 1")
  }
  if(dmg_rate > 1 || dmg_rate < 0){
    stop("dmg_rate must be a numeric between 0 and 1")
  }
  
  #determine seed value of not designated
  if(is.na(seed)){
    seed = sample.int(.Machine$integer.max, size = 1)
  }
  if(is.numeric(seed) == F){
    stop("random seed must be numeric.")
  }
  
  nsam = nrow(G)
  set.seed(seed)
  NA_seeds = sample.int(.Machine$integer.max, size = nsam)
  
  #main function operations start here
  
  #add missingness ----
  snp = ncol(G)
  n = snp*missing_rate
  for(i in 1:nsam){
    set.seed(NA_seeds[i])
    index = sample(x = snp, size = floor(n), replace = F)
    G[i,index] = NA
  }
  
  #deamination ----
  n = snp*trans_rate
  num_trans = round(n)
  
  set.seed(seed)
  #random seed to determine which columns are considered transitions
  trans_seeds = sample.int(.Machine$integer.max, size = 1) 
  #random seed to determine the coin flips for each transition
  flip_seeds = sample.int(.Machine$integer.max, size = num_trans)
  #random seed to determine which elements of the transition column actually 
  #undergo a transition. 
  deam_seeds = sample.int(.Machine$integer.max, size = snp)
  
  #randomly sample columns to become transition sites
  set.seed(trans_seeds)
  trans_index = sample(x = snp, size = num_trans, replace = F)
  
  #loop over all transition sites/columns
  for(i in 1:length(trans_index)){
    set.seed(flip_seeds[i])
    #flip coin to see if we change 0's to 1's (T), or 1's to 0's (F). 
    coin_flip = purrr::rbernoulli(1, p = 0.5)
    
    #determine elements of each column to convert
    set.seed(deam_seeds[i])
    deam_index = purrr::rbernoulli(nsam, p = dmg_rate)
    
    #find the row indices corresponding to potential elements to change
    if(coin_flip){
      possible_transitions = which(G[,i] == 0)
    } else {
      possible_transitions = which(G[,i] == 1)
    }

    #Go to next iteration if there are no possible base changes to make. 
    #Typically, this scenario should not occur. 
    if(purrr::is_empty(possible_transitions)){
      next
    }
    
    #loop over all elements in a column that could be changed.
    for(j in possible_transitions){
      if(deam_index[j] && coin_flip){
        G[j,i] = 1 
      }
      if(deam_index[j] && coin_flip==F){
        G[j,i] = 0
      }
    }
  }
  return(G)
}