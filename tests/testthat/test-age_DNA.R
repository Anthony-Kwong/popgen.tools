#Test function for age_DNA

test_age_DNA <- function(G, missing_rate, trans_rate = 0.776, dmg_rate = 0.05, seed = NA){
  nsam = nrow(G)
  snp = ncol(G)
  
  set.seed(seed)
  NA_seeds = sample.int(.Machine$integer.max, size = nsam)
  deam_seeds = sample.int(.Machine$integer.max, size = nsam)
  
  #add missingness
  n = snp*missing_rate
  for(i in 1:nsam){
    set.seed(NA_seeds[i])
    index = sample(x = snp, size = floor(n), replace = F)
    G[i,index] = NA
  }
  
  #deamination
  n = snp*trans_rate
  
  for(i in 1:nsam){
    set.seed(deam_seeds[i])
    deam_index = sample(x = snp, size = floor(n), replace = F)
    transitions = which( G[i,deam_index] == 0)
    if(length(transitions) == 0){
      next 
    }
    ntrans = length(transitions)
    trans_success = purrr::rbernoulli(ntrans, p = dmg_rate)
    
    trans_index = deam_index[transitions]
    for(k in 1:ntrans){
      if(trans_success[k]){
        G[i, trans_index[k] ] = 1
      }
    }
  }
  return(G)
}

test_that("age_DNA works",{

  #manual check with default values of trans_rate, dmg_rate
  
  set.seed(128)
  SNP=5
  missing_rate = 0.1
  seed = 1
  #G <- matrix(sample(0:1, size = SNP*2, replace = TRUE), nc = SNP)
  G <- rbind( rep(1,10),rep(0,10) )
  imp_output = age_DNA(G, missing_rate = missing_rate, seed  = seed)
  test_output = test_age_DNA(G, missing_rate = missing_rate, seed  = seed)
  
  
  nsam = nrow(G)
  set.seed(seed)
  NA_seeds = sample.int(.Machine$integer.max, size = nsam)
  deam_seeds = sample.int(.Machine$integer.max, size = nsam)
  
  #add missingness
  snp = ncol(G)
  nsam = nrow(G)
  n = snp*missing_rate
  for(i in 1:nsam){
    set.seed(NA_seeds[i])
    index = sample(x = snp, size = floor(n), replace = F)
    G[i,index] = NA
  }
  
  #add deamination
  n = snp*0.776
  for(i in 1:nsam){
    set.seed(deam_seeds[i])
    deam_index = sample(x = snp, size = floor(n), replace = F)
    valid_trans = which( G[i,deam_index] == 0)
    num_trans = length(valid_trans)
    if(num_trans == 0){
      next
    }
    trans_success = purrr::rbernoulli(num_trans, p = 0.05)
    trans_index = deam_index[valid_trans]
    for(k in 1:num_trans){
      if(trans_success[k]){
        G[i,trans_index[k]] = 1
      }
    }
  }
  
  expect_equal(G, imp_output)
  expect_equal(G, test_output)
  
  #auto check
  
  set.seed(128)
  SNP = 20 
  missing_rate = 0.2
  seed = 102
  G <- matrix(sample(0:1, size = SNP*10, replace = TRUE), nc = SNP)
  imp_output = age_DNA(G, missing_rate = missing_rate, seed = seed)
  test_output = test_age_DNA(G, missing_rate = missing_rate, seed = seed)
  expect_equal(imp_output, test_output)
  
  set.seed(128)
  SNP = 20
  missing_rate = 0.15
  seed = 1012
  G <- matrix(sample(0:1, size = SNP*20, replace = TRUE), nc = SNP)
  imp_output = age_DNA(G, missing_rate = missing_rate, trans_rate = 0.5, dmg_rate = 0.1,seed = seed)
  test_output = test_age_DNA(G, missing_rate = missing_rate, trans_rate = 0.5, dmg_rate = 0.1, seed = seed)
  expect_equal(imp_output, test_output)
})