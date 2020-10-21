test_that("sample_vec works",{
  x <- seq(1:5)
  set.seed(2)
  R_says = sample_vec(x)
  set.seed(2)
  ans = sample(x,1)
  expect_equal(ans, R_says)
  
  y <- c(1)
  R_says = sample_vec(y)
  expect_equal(1,R_says)
})