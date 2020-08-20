test_that("random_impute works",{
  #manual check
  set.seed(2)
  G = matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
  G[1,1] = NA
  G[3,4] = NA
  G[5,2] = NA
  set.seed(3)
  out = random_impute(G)
  
  G[1,1] = 0
  G[5,2] = 0
  G[3,4] = 1
  expect_equal(out, G)
})