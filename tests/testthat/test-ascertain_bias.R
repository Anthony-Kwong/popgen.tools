test_that("ascertain_bias function works",{
  #manual check
  set.seed(21)
  G = matrix(sample(0:1, size = 56, replace = TRUE), nc = 8)
  output = ascertain_bias(G, c(6,7))
  ans = G[1:5, c(1,3,5,7,8)]
  expect_equal(output, ans)
  
  set.seed(2664)
  G = matrix(sample(0:1, size = 100, replace = TRUE), nc = 10)
  output = ascertain_bias(G, c(9,10))
  ans = G[1:8, c(1,2)]
  expect_equal(output, ans)
})