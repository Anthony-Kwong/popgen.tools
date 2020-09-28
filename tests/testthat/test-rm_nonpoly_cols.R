test_that("rm_nonpoly_cols function works",{
  
  set.seed(123)
  G <-matrix(sample(0:1, size = 24, replace = TRUE), nc = 8)
  ans <- G[,-c(1,7)]
  expect_equal(ans, rm_nonpoly_cols(G)[[1]])
  expect_equal(seq(1,8)[-c(1,7)], rm_nonpoly_cols(G)[[2]])
  
  set.seed(2)
  G <- matrix(sample(0:1, size = 40, replace = TRUE), nc = 10)
  G[,1] = rep(1,4)
  G[,5] = rep(0,4)
  ans <- G[,-c(1,5)]
  expect_equal(ans, rm_nonpoly_cols(G)[[1]])
  expect_equal(seq(1,10)[-c(1,5)], rm_nonpoly_cols(G)[[2]])
  
  #check case of no polymorphic columns
  G <- rbind(c(0,0,0,0),c(0,0,0,0))
  ans <- list(NaN,NaN)
  suppressWarnings(
    expect_equal(ans,rm_nonpoly_cols(G))
  )
  
  G <- rbind(c(1,1,1,1),c(1,1,1,1))
  ans <- list(NaN,NaN)
  suppressWarnings(
    expect_equal(ans,rm_nonpoly_cols(G))
  )

  
})