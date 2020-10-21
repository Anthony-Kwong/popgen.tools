test_that("majority_flip works",{
  #results were checked manually
  set.seed(3)
  G = matrix(sample(0:1, size =20 , replace = TRUE), nc = 5)
  R_says = majority_flip(G, seed = 2)
  ans = cbind(G[,1], rep(1,4), rep(1,4), rep(1,4), G[,5])
  expect_equal(R_says, ans)
  
  set.seed(34)
  G = matrix(sample(0:1, size =35 , replace = TRUE), nc = 7)
  R_says = majority_flip(G, seed = 1)
  ans = cbind(c(0,0,0,1,0), rep(1,5), rep(1,5), 
              rep(0,5), c(0,0,0,0,1), c(0,1,0,0,0), rep(0,5))
  expect_equal(R_says, ans)
})