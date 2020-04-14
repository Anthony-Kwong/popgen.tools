test_that("have_na works",{
  M = matrix(c(1, 2, NA, 4), nrow = 2)
  expect_equal(have_na(M), T)
  
  L = matrix(NA,nrow=1)
  expect_equal(have_na(L), T)
  
  K = matrix(c(1, 0, 1, 0), nrow = 2)
  expect_equal(have_na(K), F)
})