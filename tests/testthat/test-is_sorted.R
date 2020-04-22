test_that("is_sorted works",{
  set.seed(2)
  y = runif(15, 0, 1)
  expect_equal(is_sorted(y, ascend = T), F)
  
  x = runif(10, 1, 2) %>% sort()
  expect_equal(is_sorted(x, ascend = T), T)
  
  k = runif(50, 0, 1) %>% sort(decreasing = T)
  expect_equal(is_sorted(k, ascend = F), T)
  
})