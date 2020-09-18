test_that("zero_impute works",{
  G <- rbind(c(1, 0, 1, NA, 0), c(1, NA, NA, 0, 1), c(0, 0, 1, 0, NA))
  ans <- rbind(c(1, 0, 1, 0, 0), c(1, 0, 0, 0, 1), c(0, 0, 1, 0, 0))
  expect_equal(ans,zero_impute(G))
  
  set.seed(53)
  G <-matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
  G <- age_DNA(G, missing_rate = 0.3)
  output <- zero_impute(G)
  G[is.na(G)] = 0
  expect_equal(output,G)
})