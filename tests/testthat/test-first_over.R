test_that("first_over works",{
  set.seed(19)
  x = runif (10, 0 , 1) %>% sort()
  target = 0.5
  expect_equal(first_over(x, target), 7)
  
  set.seed(1707)
  k = runif(20, 0 , 1) %>% sort()
  ans = suppressWarnings(
    first_over(k, 2)
  ) 
  expect_equal(as.numeric(NA), ans)
})