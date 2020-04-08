test_that("downsample_mat works",{
  set.seed(5)
  seq <-matrix(sample(0:1, size =40 , replace = TRUE), nc = 10)
  p <- 0.4
  seed = 52
  
  output <- downsample_mat(seq, p , seed)
  
  nSNP = ncol(seq)
  k = round (nSNP*p)
  set.seed(seed)
  kcol = sample(1:nSNP, k ) %>% sort()
  act = seq[,kcol]
  
  expect_equal(output, act)
  
  #case 2
  
  set.seed(235)
  seq <-matrix(sample(0:1, size =220 , replace = TRUE), nc = 22)
  p <- 0.7
  seed = 252
  
  output <- downsample_mat(seq, p , seed)
  
  nSNP = ncol(seq)
  k = round (nSNP*p)
  set.seed(seed)
  kcol = sample(1:nSNP, k ) %>% sort()
  act = seq[,kcol]
  
  expect_equal(output, act)
})