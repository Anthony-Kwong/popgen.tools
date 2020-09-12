test_that("het_finder works",{
  #manual checks, recall cpp uses 0 indexing
  set.seed(2)
  G = matrix(sample(0:1, size = 15, replace = TRUE), nc = 5)
  output = het_finder(G)
  expect_equal(output,c(0,3,4)+1)
  
  set.seed(2)
  G = matrix(sample(0:1, size = 40, replace = TRUE), nc = 8)
  output = het_finder(G)
  expect_equal(output,seq(0,7)+1)
  
  set.seed(3)
  G = matrix(sample(0:1, size = 40, replace = TRUE), nc = 20)
  output = het_finder(G)
  expect_equal(output,c(0,1,3,4,6,8,9,10,11,13,16,18)+1)
  
  set.seed(3)
  K = rbind(c(0,0,0),c(0,0,0))
  output = het_finder(K)
  expect_equal(output, NaN)
})