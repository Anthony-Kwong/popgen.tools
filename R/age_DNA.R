#' age_DNA function
#' 
#' A function to simulating DNA aging on discoal simulations. To simulate missingness, for each row 
#' a percentage of elements are randomly sampled and become missing. To simulate deamination
#' a percentage of elements in each row are randomly sampled according to trans_rate. Then 
#' if the element is 0, it becomes 1 with probability designated by  dmg_rate. 
#' 
#' doi: 10.1534/genetics.112.139949
#'
#' @param G : A binary, numeric genome matrix from discoal. 
#' @param missing_rate : Probability of elements that are randomly sampled to become NA for each row in G.
#' @param trans_rate: Probability of elements that are randomly sampled to be transitions for each row in G. 
#' Default is 0.776.
#' @param dmg_rate: Probability that a randomly sampled transition with an ancestral allele 0 is convert to 1.
#' Default is 0.05. 
#' @param seed: A numeric random seed. Optional.
#' @return
#' @export
#'
#' @examples
#' 
#' @importFrom purrr rbernoulli
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
  deam_seeds = sample.int(.Machine$integer.max, size = nsam)
  
  #main function operations start here
  
  #add missingness
  snp = ncol(G)
  n = snp*missing_rate
  for(i in 1:nsam){
    set.seed(NA_seeds[i])
    index = sample(x = snp, size = floor(n), replace = F)
    G[i,index] = NA
  }
  
  #deamination
  n = snp*trans_rate
  
  for(i in 1:nsam){
    #randomly sample elements for deamination
    set.seed(deam_seeds[i])
    deam_index = sample(x = snp, size = floor(n), replace = F)
    #find the ancestral alleles among sampled elements
    transitions = which(G[i,deam_index] == 0)
    #skip to next iteration if there are no value elements for transitions
    if(length(transitions) == 0){
      next 
    }
    
    #turn 0 to 1 w.p dmg_rate
    n_trans = length(transitions)
    trans_success = purrr::rbernoulli(n_trans, p = dmg_rate)
    
    #indices of elements to undergo transition
    trans_index = deam_index[transitions]
    for(j in 1:length(trans_index)){
      if(trans_success[j]){
        G[i,trans_index[j]] = 1
      }
    }
  }
  
  return(G)
}