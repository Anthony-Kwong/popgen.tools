test_that("add_missingness function works",{
  missing_rate = 0.4
  seed = 3
  
  set.seed(2)
  G = matrix(sample(0:1, size = 40, replace = TRUE), nc = 8)
  output = add_missingness(G, missing_rate = missing_rate, seed = seed)
  
  test_G = G
  set.seed(seed)
  NA_seeds = sample.int(.Machine$integer.max, size = ncol(test_G)) 
  
  for(i in 1:nrow(test_G)){
    set.seed(NA_seeds[i])
    index = purrr::rbernoulli(n = ncol(test_G), p = missing_rate)
    test_G[i,which(index == TRUE)] = NA
  }
  
  expect_equal(output, test_G)
})