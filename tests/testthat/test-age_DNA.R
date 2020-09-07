#Test function for age_DNA

test_age_DNA <- function(G, missing_rate, trans_prop = 0.776, dmg_rate = 0.05, seed = NA){
  nsam = nrow(G)
  snp = ncol(G)
  
  if(is.na(seed)){
    seed = sample.int(.Machine$integer.max, size = 1) 
  }
  
  #add missingness
  G = add_missingness(G, missing_rate = missing_rate, seed = seed)
  
  #deamination
  n_trans = snp*trans_prop
  n_trans = round(n_trans)
  
  set.seed(seed)
  #random seed to determine which columns are considered transitions
  trans_seeds = sample.int(.Machine$integer.max, size = 1) 
  #random seed to determine the coin flips for each transition
  flip_seeds = sample.int(.Machine$integer.max, size = n_trans)
  #random seed to determine which elements of the transition column actually 
  #undergo a transition. 
  deam_seeds = sample.int(.Machine$integer.max, size = snp)
  
  set.seed(trans_seeds)
  trans_sites = sample(x = snp, size = n_trans, replace =F)
  
  for (i in 1:length(trans_sites)){
    set.seed(flip_seeds[i])
    coin_flip = purrr::rbernoulli(1, p = 0.5)
    
    set.seed(deam_seeds[i])
    deam_index = purrr::rbernoulli(nsam, p = dmg_rate)
    
    if(coin_flip){
      poss_trans = which(G[,i] == 0)
    } else {
      poss_trans = which(G[,i] == 1)
    }
    
    if(purrr::is_empty(poss_trans)){
      next
    }
    
    for(j in poss_trans){
      if(deam_index[j] && coin_flip){
        G[j,i] = 1
      }
      if(deam_index[j] && coin_flip == F){
        G[j,i] = 0
      }
    }
  }
  return(G)
}

test_that("age_DNA works",{

  #manual check with default values of trans_prop, dmg_rate
  
  set.seed(128)
  SNP=10
  missing_rate = 0.1
  seed = 1
  G <- matrix(sample(0:1, size = SNP*5, replace = TRUE), nc = SNP)
  imp_output = age_DNA(G, missing_rate = missing_rate, seed  = seed)
  test_output = test_age_DNA(G, missing_rate = missing_rate, seed  = seed)
  


  #add missingness
  G = add_missingness(G, missing_rate = missing_rate, seed = seed)
  
  #add deamination
  snp = ncol(G)
  nsam = nrow(G)
  n = snp*0.776 #use default transition rate
  num_trans = round(n)
  
  set.seed(seed)
  trans_seeds = sample.int(.Machine$integer.max, size = 1) 
  flip_seeds = sample.int(.Machine$integer.max, size = num_trans)
  deam_seeds = sample.int(.Machine$integer.max, size = snp)
  
  set.seed(trans_seeds)
  trans_index = sample(x = snp, size = num_trans, replace = F)
  
  for(i in 1:length(trans_index)){
    set.seed(flip_seeds[i])
    coin_flip = purrr::rbernoulli(1, p = 0.5)
    
    set.seed(deam_seeds[i])
    deam_index = purrr::rbernoulli(nsam, p = 0.05) #use default damage rate
    
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
  
  expect_equal(G, imp_output)
  expect_equal(G, test_output)
  
  #manual check 2 ----
  
  set.seed(128)
  SNP=10
  nsam = 10
  missing_rate = 0.2
  seed = 102
  G <- matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  imp_output = age_DNA(G, missing_rate = missing_rate, seed  = seed)
  test_output = test_age_DNA(G, missing_rate = missing_rate, seed  = seed)
  
  #add missingness
  G = add_missingness(G, missing_rate = missing_rate, seed = seed)
  
  n = SNP*0.776 #use default transition rate
  num_trans = round(n)
  
  set.seed(seed)
  trans_seeds = sample.int(.Machine$integer.max, size = 1) 
  flip_seeds = sample.int(.Machine$integer.max, size = num_trans)
  deam_seeds = sample.int(.Machine$integer.max, size = SNP)
  
  set.seed(trans_seeds)
  trans_index = sample(x = SNP, size = num_trans, replace = F)
  
  for(i in 1:length(trans_index)){
    set.seed(flip_seeds[i])
    coin_flip = purrr::rbernoulli(1, p = 0.5)
    
    set.seed(deam_seeds[i])
    deam_index = purrr::rbernoulli(nsam, p = 0.05) #use default damage rate
    
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
  
  expect_equal(G, imp_output)
  expect_equal(G, test_output)
  
  #auto check ----
  
  SNP = 10 
  missing_rate = 0.2
  seed = 102
  set.seed(128)
  G <- matrix(sample(0:1, size = SNP*10, replace = TRUE), nc = SNP)
  imp_output = age_DNA(G, missing_rate = missing_rate, seed = seed)
  test_output = test_age_DNA(G, missing_rate = missing_rate, seed = seed)
  expect_equal(imp_output, test_output)
  
  set.seed(128)
  SNP = 20
  missing_rate = 0.15
  seed = 1012
  G <- matrix(sample(0:1, size = SNP*20, replace = TRUE), nc = SNP)
  imp_output = age_DNA(G, missing_rate = missing_rate, trans_prop = 0.5, dmg_rate = 0.1,seed = seed)
  test_output = test_age_DNA(G, missing_rate = missing_rate, trans_prop = 0.5, dmg_rate = 0.1, seed = seed)
  expect_equal(imp_output, test_output)
})