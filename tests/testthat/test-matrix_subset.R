test_that("matrix_subset works",{
  set.seed(2)
  G = matrix(sample(0:1, size = 15, replace = TRUE), nc = 5)
  output = matrix_subset(G,c(1,3))
  expect_equal(output, G[,c(1,3)])
  
  set.seed(4)
  G = matrix(sample(0:1, size = 60, replace = TRUE), nc = 10)
  output = matrix_subset(G,c(1))
  expect_equal(output, as.matrix(G[,c(1)]))
  
  set.seed(7)
  G = matrix(sample(0:1, size = 60, replace = TRUE), nc = 20)
  sel_cols = sort(sample.int(20,5))
  output = matrix_subset(G,sel_cols)
  expect_equal(output, G[,sel_cols])
})