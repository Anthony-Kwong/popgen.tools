test_that("any_sug works",{
  x<-c(0,1,3,5,7)
  expect_equal(any_sug(x<2),T)
  expect_equal(any_sug(x>8),F)
})