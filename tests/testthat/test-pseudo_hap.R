test_that("pseudo_hap works",{
  
  #manual check
  set.seed(122)
  SNP=5
  missing_rate = 0.1
  seed = 12311
  G <- matrix(sample(0:1, size = SNP*4, replace = TRUE), nc = SNP)
  output = pseudo_hap(G, seed)
  
  nsam = nrow(G)
  snp = ncol(G)
  set.seed(seed)
  p_seeds = sample.int(.Machine$integer.max, size = nsam/2)
  
  start_indices = seq(1, nsam, by = 2)
  new_haps = list()
  for(i in 1:length(start_indices)){
    set.seed(p_seeds[i])
    choice = sample(c(0,1) , replace=TRUE, size= snp)
    new_row = rep(NA, snp)
    for (j in 1:snp){
      if (choice [j] ==1){
        new_row[j] = G[start_indices[i], j] 
      } else {
        new_row[j] = G[start_indices[i], j] 
      }
    }
    new_haps[[i]] = new_row
  }
  actual = rlist::list.rbind(new_haps)
  
  expect_equal(actual, output)
})